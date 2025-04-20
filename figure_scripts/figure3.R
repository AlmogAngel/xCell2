library(tidyverse)
library(survival)
library(patchwork)
library(survminer)
library(shapviz)

source("~/xCell2_dev/paper_figures/load_lightgbm_functions.R")

set.seed(123)


# Load data
icb_expr <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_expression_tpm.rds")
metadata <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_metadata.rds")
all_methods_predictions <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_all_methods_predictions.rds")


# ----- Prepare metadata ----

categorical_features <- c('treatment', 'Cancer_Type', 'Gender', 'Response')
metadata <- metadata[,categorical_features]
metadata <- na.omit(metadata)
metadata$treatment <- as.factor(metadata$treatment)
metadata$Cancer_Type <- as.factor(metadata$Cancer_Type)
metadata$Gender <- as.factor(metadata$Gender)

# Ensure Response is numeric feature (0/1) for LightGBM.
metadata$Response <- ifelse(metadata$Response == "Response", 1, 0)
categorical_features <- categorical_features[categorical_features != 'Response']

# ----- Prepare training data ----

method2use <- c("xcell2", "bayesprism", "cbrx", "decon", "epic", "mcpcounter", "dtangle", "scaden", "dwls", "bisque", "scdc")
refs2use <- c("tme_c", "sc_pan_cancer")

methods_data <- list()
for (m in method2use) {

  ref1_data <- all_methods_predictions[[m]][[refs2use[1]]] %>% t() %>% as.data.frame()
  ref2_data <- all_methods_predictions[[m]][[refs2use[2]]] %>% t() %>% as.data.frame()
  
  if (all(rownames(ref1_data) == rownames(ref2_data))) {
    
    colnames(ref1_data) <- paste0(colnames(ref1_data), "_", refs2use[1])
    colnames(ref2_data) <- paste0(colnames(ref2_data), "_", refs2use[2])
    
    data <- cbind(ref1_data, ref2_data)
    colnames(data) <- make.unique(colnames(data))
    
    colnames(data) <- colnames(data) %>%
      gsub("^X", "", .) %>%
      gsub("\\.|\\-|\\s|\\+", "_", .) %>%
      gsub("#", "_", .) %>%
      gsub(",", "_", .)

    methods_data[[m]] <- data
    
  }else{
    errorCondition("Samples doesn't matach")
  }
  
}

# Round xCell2 results
methods_data$xcell2 <- round(methods_data$xcell2, 3)

# Load TIDE
tide = read.csv("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/dataTIDE.csv", header = TRUE) %>% as.data.frame()
rownames(tide) <- tide[,1]
tide <- tide[,-1]
methods_data[["tide"]] <- tide[rownames(methods_data[[1]]),]

# Load EASIER
easier <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_easier.rds") %>% t() %>% as.data.frame()
easier <- t(easier)
methods_data[["easier"]] <- easier[rownames(methods_data[[1]]),]

# Load IMPERS
impers <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_impres.rds")
a <- as.data.frame(impers[rownames(methods_data[[1]]),])
rownames(a) <- rownames(impers)
methods_data[["impers"]] <- a

# Add clinical data (metadata) as features for all methods
methods_data <- lapply(methods_data, function(m){
  m <- m[rownames(metadata),]
  cbind(m, metadata)
})

methods_data[["clinical_data"]] <- metadata

method2use <- c(method2use, "tide", "easier", "impers", "clinical_data")

# ------ Run LightGBM model -----

for (method in method2use) {
  
  outfile <- paste0("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/lightgbm_res_", method, ".rds")
  if (!file.exists(outfile)) {
    print(paste0("Running LightGBM for: ", method, "..."))
    res <- run_lightgbm(data = methods_data[[method]], outer_iterations = 100, ncores = 60)
    saveRDS(res, outfile)
  }
  
}


f <- list.files(path = "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/", pattern = "lightgbm_res*", full.names = TRUE)
ff <- lapply(f, function(i){
  i_in <- readRDS(i)
  i_in$auc_values
})
names(ff) <- basename(f)
names(ff) <- gsub(".rds", "", gsub("lightgbm_res_", "", names(ff)))
sort(sapply(ff, mean), decreasing = TRUE)
sort(sapply(ff, median), decreasing = TRUE)


# Plot ---------


ff <- bind_rows(ff) |> 
  pivot_longer(cols = everything(), names_to = "method", values_to = "auc")

# Ensure method is a factor with specified levels
ff <- ff %>%
  # mutate(method = factor(method, levels = c("xcell2", "cbrx", "dtangle",  "tide", "decon", "easier",
  #                                           "bayesprism", "mcpcounter", "epic", "impers", "metadata")),
  #        is_xcell2 = ifelse(method == "xcell2", "yes", "no"))
mutate(is_xcell2 = ifelse(method == "xcell2", "yes", "no"))

# Perform statistical tests for each grid group
xcell2_vs_cbrx <- ggpubr::compare_means(auc ~ method, data = grid_res %>% filter(method %in% c("xcell2", "mcpcounter")),
                                        method = "wilcox.test",
                                        group.by = "grid")
xcell2_vs_metadata <- ggpubr::compare_means(auc ~ method, data = grid_res %>% filter(method %in% c("xcell2", "clinical_data")),
                                            method = "wilcox.test",
                                            group.by = "grid")


ggplot(ff, aes(x = reorder(method, -auc, FUN = median), y = auc, fill = is_xcell2)) +
  geom_boxplot(width = .5, show.legend = F, position = "dodge", outlier.colour=NA) +
  scale_fill_manual(values = c("yes" = "tomato", "no" = "gray")) +
  theme_minimal() +
  labs(y = "AUC", x="") +
  scale_y_continuous(limits = c(0.55, 0.72), breaks = seq(0.5, 0.73, by = 0.05)) +
  scale_x_discrete(labels = c(
    "xcell2" = "xCell2",
    "dtangle" = "dtangle",
    "cbrx" = "CIBERSORTx",
    "bayesprism" = "BayesPRISM",
    "decon" = "DeconRNASeq",
    "epic" = "EPIC",
    "mcpcounter" = "MCPcounter",
    "clinical_data" = "Clinical data",
    "tide" = "TIDE",
    "impers" = "IMPERS",
    "easier" = "EASIER"
  )) +
  ggpubr::stat_compare_means(comparisons = list(c("xcell2", "mcpcounter"), c("xcell2", "clinical_data")),
                             method = "wilcox.test",
                             label = "p.format",
                             label.y = c(0.68, 0.7)) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1), 
    axis.title.y = element_text(size = 20)
  )


# Survival analysis ------------

surv_plot <- function(data, strat, plot_title) {
  
  dat <- data[!is.na(data[[strat]]),]
  #dat$time <- dat$time/365
  labs <- sort(unique(dat[[strat]]))
  
  labels2use <- unique(dat[,strat]) %>% pull()
  xlim2use <- min(max(dat[dat[[strat]] == labels2use[1],]$time),
                  max(dat[dat[[strat]] == labels2use[2],]$time))
  xlim2use <- floor(xlim2use / 5) * 5
  
  splot <- survminer::ggsurvplot(fit = surv_fit(as.formula(paste0("Surv(time, status) ~", strat)), data = dat),
                                 data = dat,
                                 palette = c("#2E8B57", "tomato"),
                                 ggtheme = theme_bw(),
                                 xlim = c(0, xlim2use),
                                 break.time.by = 5,
                                 pval = TRUE,
                                 pval.method = TRUE,
                                 pval.size = 3.5,
                                 pval.coord = c(0,0.05),
                                 pval.method.coord = c(0,0.12),
                                 risk.table = TRUE,
                                 risk.table.y.text = FALSE,
                                 risk.table.height = 0.3,
                                 tables.y.text = TRUE,
                                 fontsize = 3,
                                 size = 0.2,
                                 conf.int = TRUE,
                                 conf.int.alpha = 0.2,
                                 conf.int.style = "ribbon",
                                 censor = T,
                                 censor.shape = 124,
                                 censor.size = 2,
                                 legend.labs = labs)
  
  splot$plot <- splot$plot +
    labs(title = paste(plot_title), fill = NULL, color = NULL,
         x = "Time (months)") +
    scale_y_continuous(labels = scales::percent_format())+
    guides(color = guide_legend(override.aes = list(shape = NA,
                                                    linewidth=0.6)))+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color="black", fill=NA,
                                      linewidth=0.6),
          legend.justification = c(1,1),
          legend.key.size = unit(4,"mm"),
          legend.spacing.y = unit(0,"mm"),
          legend.background = element_rect(color="black", fill=NA,
                                           linewidth = 0.6*0.5),
          legend.position = c(1,1),
          legend.box.margin = margin(t = 0.6*(.pt*72.27/96/4),
                                     r = 0.6*(.pt*72.27/96/4)),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color="black"),
          axis.title = element_text(size = 16, color="black"),
          plot.background = element_rect(fill="white"),
          plot.tag.position = c(0,1))
  
  
  splot$table <- splot$table +
    labs(title="No. at risk", x = NULL, y = NULL) +
    coord_cartesian(clip = "off", xlim=c(0,xlim2use),ylim=c(0,2)) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size=10),
          panel.border = element_rect(color = NA),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  splot
}


surv.in <- readRDS("/bigdata/loainaom/overallSurvival.rds")
metadata <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_metadata.rds")
surv.in$Sample <- rownames(surv.in)
metadata$Sample <- rownames(metadata)


f <- list.files(path = "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/", pattern = "lightgbm_res*", full.names = TRUE)
all_surv_data <- lapply(f, function(i){
  i_in <- readRDS(i)
  i_in <- lapply(i_in$model_predictions, function(z){z[1]$predictions})
  i_in <- bind_rows(i_in)
  
  surv_data <- i_in %>%
    pivot_longer(cols = everything(), names_to = "Sample", values_to = "Prediction") %>%
    group_by(Sample) %>%
    summarise(mean_pred = mean(Prediction, na.rm = TRUE)) %>%
    mutate(xcell2_classification = ifelse(mean_pred <= median(mean_pred), "Below median", "Above median"))
  
  surv_data <- surv_data %>%
    right_join(surv.in[,c("survival_time", "survival_status", "Sample")], by = "Sample") %>%
    right_join(metadata[,c("Cancer_Type", "Gender", "treatment", "Sample")])
  
  colnames(surv_data)[4] <- "time"
  colnames(surv_data)[5] <- "status"
  surv_data$Cancer_Type <- as.factor(surv_data$Cancer_Type)
  
  surv_data$status <- surv_data$status+1
  
  surv_data
})
names(all_surv_data) <- basename(f)
names(all_surv_data) <- gsub(".rds", "", gsub("lightgbm_res_", "", names(all_surv_data)))

surv_data.list.xcell2 <- group_split(all_surv_data$xcell2, Cancer_Type) %>%
  setNames(surv_data %>% group_keys(Cancer_Type) %>%
             pull(Cancer_Type))

pl <- list(surv_plot(data = surv_data.list.xcell2$Melanoma, strat = "xcell2_classification",
                     plot_title = "Melanoma"),
           surv_plot(data = surv_data.list.xcell2$NSCLC, strat = "xcell2_classification",
                     plot_title = "NSCLC"),
           surv_plot(data = surv_data.list.xcell2$`Urothelial Carcinoma`, strat = "xcell2_classification",
                     plot_title = "Urothelial Carcinoma"))

pl[[2]]$plot <- pl[[2]]$plot + labs(y = NULL)
pl[[3]]$plot <- pl[[3]]$plot + labs(y = NULL)
pl[[1]]$plot <- pl[[1]]$plot + labs(x = NULL)
pl[[3]]$plot <- pl[[3]]$plot + labs(x = NULL)


p <- pl[[1]]$plot + pl[[2]]$plot + pl[[3]]$plot +
  #  pl[[1]]$table + pl[[2]]$table + pl[[3]]$table +
  plot_layout(design = "ABC\nDEF", height=c(1,0.3))

# Cox Proportional Hazards Model ------------


surv_data.xcell2 <- all_surv_data$xcell2 |> 
  na.omit()

surv_obj <- Surv(time = surv_data.xcell2$time, event = surv_data.xcell2$status)

cox_fit <- coxph(surv_obj ~ mean_pred + Cancer_Type + Gender + treatment, data = surv_data.xcell2)
summary(cox_fit)

#  SHAP ------------------


res <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/lightgbm_res_xcell2.rds")
res <- lapply(res$shap_values_list, function(x){as.data.frame(get_shap_values(x))})
res <- as_tibble(bind_rows(res, .id = "model_id"), rownames = "sample")
res$sample <- gsub("\\...*.", "", res$sample)

importance <- res %>%
  pivot_longer(cols = -c(sample, model_id), names_to = "feature", values_to = "shap") %>%
  group_by(feature, sample) |> 
  summarise(MeanAbsSHAP = abs(mean(shap)))

top10_features <- importance |> 
  group_by(feature) |> 
  summarise(mediaMeanAbsSHAP = median(MeanAbsSHAP)) |> 
  arrange(desc(mediaMeanAbsSHAP)) %>%
  slice_head(n = 10) |> 
  pull(feature)

importance_top10 <- importance |> 
  dplyr::filter(feature %in% top10_features)
  
importance_top10$ref <- NA
importance_top10$ref[grepl("sc_pan_cancer", importance_top10$feature)] <- "Pan Cancer"
importance_top10$ref[grepl("tme_c", importance_top10$feature)] <- "TME Compendium"

importance_top10$feature[importance_top10$feature == "B_cell_sc_pan_cancer"] <- "B cell"
importance_top10$feature[importance_top10$feature == "CD8_positive__alpha_beta_T_cell_tme_c"] <- "CD8+ T cell"
importance_top10$feature[importance_top10$feature == "T_helper_1_cell_tme_c"] <- "T helper 1"
importance_top10$feature[importance_top10$feature == "T_helper_2_cell_tme_c"] <- "T helper 2"
importance_top10$feature[importance_top10$feature == "CD8__T_cell_PD1_high_tme_c"] <- "CD8+ T cell (PD1 high)"
importance_top10$feature[importance_top10$feature == "fibroblast_tme_c"] <- "Fibroblast"
importance_top10$feature[importance_top10$feature == "inflammatory_macrophage_tme_c"] <- "Inflammatory Macrophage"
importance_top10$feature[importance_top10$feature == "malignant_cell_tme_c"] <- "Malignant cell"
importance_top10$feature[importance_top10$feature == "naive_thymus_derived_CD4_positive__alpha_beta_T_cell_sc_pan_cancer"] <- "CD4+ T cell (naive)"
importance_top10$feature[importance_top10$feature == "myeloid_cell_tme_c"] <- "Myeloid cell"

importance_top10$feature <- as.factor(importance_top10$feature)
importance_top10$ref <- as.factor(importance_top10$ref)

color_palette <- c(
  "#FBB4AE", "#33ccd0"
)


ggplot(importance_top10, aes(x = reorder(feature, MeanAbsSHAP, FUN = median), y = MeanAbsSHAP, fill = ref)) +
  geom_boxplot(width = 0.6, alpha=0.9) +
  scale_fill_manual(values = color_palette) +
  coord_flip(ylim = c(0, 0.25)) +
  labs(x = "",
       y = "mean(|SHAP|)",
       fill = "Reference:") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 18,margin = margin(15)),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title  = element_text(size = 18, face = "bold"),
    legend.position = "top"
  )


