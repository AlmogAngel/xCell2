library(tidyverse)
source("/bigdata/almogangel/xCell2/dev_scripts/loai/functions_loai_almog.R")
set.seed(123)


# ------ load data and metadata -----
TPMs <- readRDS("/bigdata/almogangel/xCell2/dev_scripts/loai/TPMs.rds")
metadata <- readRDS("/bigdata/almogangel/xCell2/dev_scripts/loai/totaldata_mini.rds")
training_data_dir <- "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/"

# ------ Run xCell2 -----

if (FALSE) {
  # kass_blood
  ref.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/kass_tumor_ref.rds")
  xcell2_object <- xCell2::xCell2Train(mix = TPMs, ref = ref.in$ref, labels = ref.in$labels, lineage_file = ref.in$lineage_file, num_threads = 40, ref_type = "rnaseq",
                                       use_sillover = TRUE, return_analysis = FALSE)
  xcell2_res <- xCell2::xCell2Analysis(mix = TPMs, xcell2object = xcell2_object, spillover = TRUE, spillover_alpha = 0.25, raw_scores = FALSE, num_threads = 40, ref_is_sc = FALSE)
  saveRDS(xcell2_res, paste0(training_data_dir, "icb_kass_tumor.xcell2.rds"))

  # sc_pan_cancer
  ref.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sc_pan_cancer_ref.rds")
  xcell2_object <- xCell2::xCell2Train(mix = TPMs, ref = ref.in$ref, labels = ref.in$labels, lineage_file = ref.in$lineage_file, num_threads = 40, ref_type = "sc",
                                       use_sillover = TRUE, return_analysis = FALSE)
  xcell2_res <- xCell2::xCell2Analysis(mix = TPMs, xcell2object = xcell2_object, spillover = TRUE, spillover_alpha = 0.25, raw_scores = TRUE, num_threads = 40, ref_is_sc = TRUE)
  saveRDS(xcell2_res, paste0(training_data_dir, "icb_sc_pan_cancer.xcell2.rds"))
}


# ----- Set run parameters -------

n_models <- 100
n_cores <- 10
add_genes <- FALSE
return_shap <- FALSE
round_xcell2_res <- 3

# label2use <- c("PD", "CR", "Response") # labels: For RECIST: "PD" "SD" "PR" "CR" or "NoResponse"/"Response"
cancers <- c("Melanoma", "Renal Cell Carcinoma", "Urothelial Carcinoma", "NSCLC")
# cancers2use <- unlist(lapply(1:length(cancers), function(x) combn(cancers, x, simplify = FALSE)), recursive = FALSE)
refs <- c("kass_tumor", "sc_pan_cancer")
# refs2use <- unlist(lapply(1:length(refs), function(x) combn(refs, x, simplify = FALSE)), recursive = FALSE)
# remove_sd <- c(TRUE, FALSE)
# params_grid <- expand.grid(label2use = label2use, cancers2use = cancers2use, refs2use = refs2use, remove_sd = remove_sd)


params_grid <- expand.grid(label2use = "PD", cancers2use = list(cancers), refs2use = list(refs), remove_sd = FALSE)


grid_res <- lapply(1:nrow(params_grid), function(i){


  # ----- Prepare metadata ----
  metadata <- dummy_metadata(metadata, cancers = params_grid$cancers2use[i][[1]], label = params_grid$label2use[i], remove_sd = params_grid$remove_sd[i])

  # ----- Prepare training datasets ----
  # Get data
  method2use <- c("xcell2", "bayesprism", "cbrx", "decon", "epic", "mcpcounter", "dtangle")

  ref2load <- params_grid$refs2use[i][[1]]
  methods_data <- list()
  for (m in method2use) {
    print(m)
    if (length(ref2load) == 1) {
      file <- list.files(training_data_dir, full.names = TRUE, pattern = paste0("icb_", params_grid$refs2use[i][[1]][1], ".", m))
      ref_data <- readRDS(file) %>% t() %>% as.data.frame()
      colnames(ref_data) <- clean_names(colnames(ref_data))
      methods_data[[m]] <- ref_data

    }else{

      file1 <- list.files(training_data_dir, full.names = TRUE, pattern = paste0("icb_", params_grid$refs2use[i][[1]][1], ".", m))
      file2 <- list.files(training_data_dir, full.names = TRUE, pattern = paste0("icb_", params_grid$refs2use[i][[1]][2], ".", m))

      ref1_data <- readRDS(file1) %>% t() %>% as.data.frame()
      ref2_data <- readRDS(file2) %>% t() %>% as.data.frame()

      if (all(rownames(ref1_data) == rownames(ref2_data))) {

        data <- cbind(ref1_data, ref2_data)
        colnames(data) <- make.unique(colnames(data))
        colnames(data) <- clean_names(colnames(data))

        methods_data[[m]] <- data

      }else{
        errorCondition("rownames doen't matach")
      }
    }

  }

  # Round xCell2 results?
  if (!is.null(round_xcell2_res)) {
    methods_data$xcell2 <- round(methods_data$xcell2, round_xcell2_res)
  }

  tide = read.csv('/bigdata/loainaom/newRuns_0708/dataFILES/dataTIDE.csv', header = TRUE) %>% as.data.frame()
  rownames(tide) <- tide[,1]
  tide <- tide[,-1]
  methods_data[["tide"]] <- tide[rownames(methods_data[[1]]),]
  easier <- readRDS('/bigdata/almogangel/xCell2_data/ICB_prediction/icb_easier.rds') %>% t() %>% as.data.frame()
  easier <- t(easier)
  methods_data[["easier"]] <- easier[rownames(methods_data[[1]]),]
  impers <- readRDS('/bigdata/almogangel/xCell2_data/ICB_prediction/icb_impres.rds')
  a <- as.data.frame(impers[rownames(methods_data[[1]]),])
  rownames(a) <- rownames(impers)
  methods_data[["impers"]] <- a

  if (add_genes) {
    genes = readxl::read_xlsx('/bigdata/almogangel/xCell2/dev_scripts/loai/genes_to_use_for_ml.xlsx') # Genes as features
    wanted_genes = genes$genes
    genes = t(TPMs[wanted_genes, rownames(metadata)]) %>% as.data.frame()
    colnames(genes) = paste0("gene_", colnames(genes))
    metadata <- cbind(genes, metadata)
  }


  # Combine metadata with methods data
  methods_data <- lapply(methods_data, function(m){
    m <- m[rownames(metadata),]
    cbind(m, metadata)
  })

  methods_data[["metadata"]] <- metadata

  method2use <- c(method2use, "tide", "easier", "impers", "metadata")

  # ------ Run LightGBM model -----

  if (return_shap) {
    method2use <- "xcell2"

    methods_shap.res <- parallel::mclapply(1:n_models, function(i){

      shap_data <- predict_response_lightgbm(data = methods_data[["xcell2"]],
                                             num_threads = 1,
                                             return_shap = return_shap)

    }, mc.cores = n_cores)

    shap_data_comb <- methods_shap.res[[1]]
    shap_data_comb$X <- t(as.data.frame(apply(shap_data_comb$X, 2, mean)))
    shap_data_comb$S <- t(as.data.frame(apply(shap_data_comb$S, 2, mean)))


    for (i in 2:length(methods_shap.res)) {
      shap_data2 <- methods_shap.res[[i]]
      cols <- intersect(colnames(shap_data_comb$X), colnames(shap_data2$X))
      X2 <- apply(shap_data2$X[,cols], 2, mean)
      shap_data_comb$X <- rbind(shap_data_comb$X[,cols], X2)

      cols <- intersect(colnames(shap_data_comb$S), colnames(shap_data2$S))
      S2 <- apply(shap_data2$S[,cols], 2, mean)
      shap_data_comb$S <- rbind(shap_data_comb$S[,cols], S2)
      shap_data_comb$baseline <- c(shap_data_comb$baseline, shap_data2$baseline)
    }
    shap_data_comb$baseline <- mean(shap_data_comb$baseline)



    shapviz::sv_importance(shap_data_comb, kind = "beeswarm")


  }else{

    tictoc::tic(paste0("Training ", n_models, " models repreats..."))
    methods_auc.res <- parallel::mclapply(1:n_models, function(i){

      methods_auc <- sapply(method2use, function(m){
        data <- methods_data[[m]]
        auc <- predict_response_lightgbm(data = data,
                                         num_threads = 1, # Do not change
                                         return_shap = return_shap)
      })
      methods_auc

    }, mc.cores = n_cores)
    tictoc::toc()

    methods_auc.res <-  bind_rows(methods_auc.res)
    return(methods_auc.res)
  }

})

saveRDS(grid_res, "/bigdata/almogangel/xCell2_data/ICB_prediction/prediction_grid_res.rds")


grid_res <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/prediction_grid_res.rds")
grid_res <- grid_res[[1]]
grid_res = grid_res %>% pivot_longer(cols = everything(), names_to = "method", values_to = "auc")
grid_res <- bind_rows(grid_res, .id = "grid")


grid_res %>%
  # filter(grid %in% 94) %>%
  ggplot(., aes(x=method, y=auc, fill=method)) +
  geom_boxplot()


# Ensure method is a factor with specified levels
grid_res <- grid_res %>%
  mutate(method = factor(method, levels = c("xcell2", "cbrx", "dtangle",  "tide", "decon", "easier",
                                            "bayesprism", "mcpcounter", "epic", "impers", "metadata")),
         is_xcell2 = ifelse(method == "xcell2", "yes", "no"))

# Perform statistical tests for each grid group
xcell2_vs_cbrx <- ggpubr::compare_means(auc ~ method, data = grid_res %>% filter(method %in% c("xcell2", "cbrx")),
                                   method = "wilcox.test",
                                   group.by = "grid")
xcell2_vs_metadata <- ggpubr::compare_means(auc ~ method, data = grid_res %>% filter(method %in% c("xcell2", "metadata")),
                                    method = "wilcox.test",
                                    group.by = "grid")

# Plot the data with p-values
ggplot(grid_res, aes(x = method, y = auc, fill = is_xcell2)) +
  geom_boxplot(width = .5, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes" = "tomato", "no" = "gray")) +
  theme_minimal() +
  labs(y = "AUC", x="") +
  scale_y_continuous(limits = c(0.58, 0.80), breaks = seq(0.6, 0.80, by = 0.05)) +
  scale_x_discrete(labels = c(
    "xcell2" = "xCell2",
    "dtangle" = "dtangle",
    "cbrx" = "CIBERSORTx",
    "bayesprism" = "BayesPRISM",
    "decon" = "DeconRNASeq",
    "epic" = "EPIC",
    "mcpcounter" = "MCPcounter",
    "metadata" = "Metadata",
    "tide" = "TIDE",
    "impers" = "IMPERS",
    "easier" = "EASIER"
  )) +
  ggpubr::stat_compare_means(comparisons = list(c("xcell2", "cbrx"), c("xcell2", "metadata")),
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = c(0.76, 0.78)) +
  theme(
    axis.text.x = element_text(size = 14), # Increase x-axis label size
    axis.title.y = element_text(size = 20)  # Increase y-axis label size
  )





surv.in <- readRDS("/bigdata/loainaom/overallSurvival.rds")

method2use <- c("xcell2", "bayesprism", "decon", "epic", "mcpcounter", "dtangle")
ref2load <- "kass_tumor"
methods_data <- list()
for (m in method2use) {
  file <- list.files(training_data_dir, full.names = TRUE, pattern = paste0("icb_", ref2load, ".", m))
  ref_data <- readRDS(file) %>% t() %>% as.data.frame()
  colnames(ref_data) <- clean_names(colnames(ref_data))
  methods_data[[m]] <- ref_data
}

method_scores <- methods_data[[1]]
method_scores <- method_scores[rownames(surv.in), "CD8__T_cell_PD1_high"]
names(method_scores) <- rownames(surv.in)

q25 <- quantile(method_scores, 0.25)
q75 <- quantile(method_scores, 0.75)
qBottom <- names(method_scores)[method_scores <= q25]
qTop <- names(method_scores)[method_scores >= q75]

surv_data <- surv.in[c(qBottom, qTop), 1:2]
surv_data$is_top <- NA
surv_data[qBottom,]$is_top <- "low"
surv_data[qTop,]$is_top <- "high"

library(survival)
library(survminer)

surv_data$survival_time <- as.numeric(surv_data$survival_time)
surv_data$survival_status <- as.numeric(surv_data$survival_status)
surv_data$is_top <- as.factor(surv_data$is_top)

surv_obj <- Surv(time = surv_data$survival_time, event = surv_data$survival_status)

# Fit the Kaplan-Meier model
fit <- survfit(surv_obj ~ is_top, data = surv_data)


ggsurvplot(fit, data = surv_data, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, legend.labs = c("Low", "Top"),
           xlab = "Time (days)", ylab = "Survival probability",
           ggtheme = theme_minimal())
