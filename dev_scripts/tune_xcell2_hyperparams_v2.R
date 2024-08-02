library(tidyverse)

set.seed(123)

methods2use <- c("xCell2", "BayesPrism", "CIBERSORTx", "DeconRNASeq", "EPIC", "MCPcounter", "dtangle")

refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))
refsRDSList <- lapply(refList, function(ref_type){
  refs <- lapply(ref_type, function(ref){
    # Load reference
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
    ref.in
  })
  names(refs) <- ref_type
  refs
})

vals.refs.res <- refval.tbl %>%
  mutate(n_val_samples = ncol(cyto.vals$truth[[val_type]][[val_dataset]])) %>%
  filter(n_shared_celltypes > 2) %>%
  mutate(method = "xCell2", .before = everything())

files <- list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations", pattern = ".cyto.res.rds", full.names = TRUE)
cyto.Res <- lapply(files, function(f){
  readRDS(f) %>%
    dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
    group_by(across(-ncol(.))) %>%
    summarise(res = list(do.call(rbind, res)), .groups = 'drop')
}) %>%
  do.call(rbind, .)

cyto.Res <- cyto.Res %>%
  rowwise() %>%
  mutate(res = list(round(res, 3)))

ref_val_pairs <- cyto.Res %>%
  group_by(ref_tissue, ref_name, val_type, val_dataset) %>%
  summarise(n = n()) %>%
  filter(n == length(methods2use)-1) %>%
  dplyr::select(-n) %>%
  ungroup()

# Random ?
# train_ds <- unique(ref_val_pairs$val_dataset)
# train_ds <- train_ds[sample(length(train_ds), length(train_ds)*0.25)]


# Choose 6 datasets that xCell 2.0 preforms worst (remove) -----------------

getCors <- function(ref, res, truth, shared_cts, shared_samp, get_spill_cors){

  celltypes <- intersect(rownames(res), rownames(truth))
  celltypes <- celltypes[celltypes %in% shared_cts]

  samples <- intersect(colnames(res), colnames(truth))
  samples <- samples[samples %in% shared_samp]


  df <- lapply(celltypes, function(ct){

    truth <- truth[ct,samples]
    res <- res[ct,samples]


    tibble(celltype = ct, truth = truth, prediction = res)

  }) %>%
    bind_rows()


  if (get_spill_cors) {
    # Generate all pairs of cell types
    celltype_pairs <- crossing(predicted_celltype = unique(df$celltype), true_celltype = unique(df$celltype))

    #     source("/bigdata/almogangel/xCell2/dev_scripts/xCell2Train_tmp_version.R")
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
    dep_list <- getDependencies(ref.in$lineage_file)

    celltype_pairs <- celltype_pairs %>%
      rowwise() %>%
      mutate(true_celltype = ifelse(true_celltype %in% c(predicted_celltype, unlist(dep_list[predicted_celltype])), NA, true_celltype)) %>%
      drop_na()

    cors.out <- celltype_pairs %>%
      rowwise() %>%
      mutate(cor = cor(
        df %>% filter(celltype == predicted_celltype) %>% pull(prediction),
        df %>% filter(celltype == true_celltype) %>% pull(truth),
        method = "spearman", use = "pairwise.complete.obs"),
        p_value = cor.test(
          df %>% filter(celltype == predicted_celltype) %>% pull(prediction),
          df %>% filter(celltype == true_celltype) %>% pull(truth),
          method = "spearman", exact = FALSE)$p.value,
        n = sum(!is.na(df %>% filter(celltype == predicted_celltype) %>% pull(prediction)) &
                  !is.na(df %>% filter(celltype == true_celltype) %>% pull(truth)))) %>%
      mutate(cor = ifelse(is.na(cor), 0, cor)) %>%
      ungroup()

    cors.out %>%
      mutate(celltype = paste0(predicted_celltype, "_", true_celltype)) %>%
      dplyr::select(celltype, cor, p_value, n) %>%
      return(.)

  }else{
    df %>%
      group_by(celltype) %>%
      dplyr::summarize(
        cor = cor(truth, prediction, method = "spearman", use = "pairwise.complete.obs"),
        p_value = cor.test(truth, prediction, method = "spearman", exact = FALSE)$p.value,
        n = sum(!is.na(truth) & !is.na(prediction))
      ) %>%
      mutate(cor = ifelse(is.na(cor), 0, cor)) %>%
      return(.)
  }

}

combineRhos <- function(rhos, sample_sizes = NULL, use_median = TRUE, summarize = FALSE){

  if (length(rhos) == 1) {
    return(rhos)
  }

  if (is.null(sample_sizes)) {
    sample_sizes <- rep(1, length(rhos))
  }


  if (length(sample_sizes) != length(rhos)) {
    # sample_sizes <- rep(sample_sizes, length(rhos))
    stop("values and weights must have the same length")
  }


  rhos[rhos == 1] <- 0.999999999
  rhos[rhos == -1] <- -0.999999999

  # Fisher's Z Transformation
  z_values <- 0.5 * log((1 + rhos) / (1 - rhos))

  if (!summarize) {
    z_weighted <- sample_sizes * z_values / mean(sample_sizes)
    rho_weighted <- (exp(2 * z_weighted) - 1) / (exp(2 * z_weighted) + 1)

    return(rho_weighted)
  }

  if (use_median) {
    weighted_median <- function(values, weights) {
      if (length(values) != length(weights)) {
        stop("values and weights must have the same length")
      }

      # Sort values and weights by values
      order_index <- order(values)
      values <- values[order_index]
      weights <- weights[order_index]

      # Calculate the cumulative sum of weights
      cum_weights <- cumsum(weights)
      total_weight <- sum(weights)

      # Find the index where the cumulative weight exceeds half of the total weight
      median_index <- which(cum_weights >= total_weight / 2)[1]

      # Return the corresponding value
      return(values[median_index])
    }
    z_median <- weighted_median(z_values, sample_sizes)
    rho_weighted_median <- (exp(2 * z_median) - 1) / (exp(2 * z_median) + 1)

    return(rho_weighted_median)
  }


  # Weighted Average of Z values
  weights <- sample_sizes
  z_mean <- sum(weights * z_values) / sum(weights)

  # Back Transformation
  rho_weighted_mean <- (exp(2 * z_mean) - 1) / (exp(2 * z_mean) + 1)
  return(rho_weighted_mean)

}

bold_xCell2_labels <- function(labels) {
  lapply(labels, function(label) {
    modified_label <- str_remove(label, "#.*")  # Remove everything after "#"
    if (startsWith(modified_label, "xCell2")) {
      bquote(bold(.(modified_label)))
    } else {
      modified_label
    }
  })
}

methods2use <- c("xCell2", "BayesPrism", "CIBERSORTx", "DeconRNASeq", "EPIC", "MCPcounter", "dtangle")



y <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/xcell2_analysis_benchmark_output.rds")
vals.refs.res$res <- y


all_cors <- parallel::mclapply(1:nrow(ref_val_pairs), function(i){


  valType <- ref_val_pairs[i,]$val_type
  valDataset <- ref_val_pairs[i,]$val_dataset
  refName <- ref_val_pairs[i,]$ref_name
  truth_mat <- cyto.vals$truth[[valType]][[valDataset]]


  cyto.Res.tmp <- cyto.Res %>%
    filter(ref_name == refName & val_dataset == valDataset)

  yy <- cyto.Res.tmp[cyto.Res.tmp$method == "BayesPrism",]
  yy$method <- "xCell2"
  yy$res <- pull(filter(vals.refs.res, ref_name == refName & val_dataset == valDataset), res)
  cyto.Res.tmp <- rbind(yy, cyto.Res.tmp)

  shared_celltypes <- Reduce(intersect, lapply(cyto.Res.tmp$res, rownames))
  shared_samples <- Reduce(intersect, lapply(cyto.Res.tmp$res, colnames))

  if (length(shared_celltypes) != unique(cyto.Res.tmp$n_shared_celltypes)) {
    errorCondition("Some method is missing cell type(s)")
  }

  out <- cyto.Res.tmp %>%
    rowwise() %>%
    mutate(cors = list(getCors(ref = ref_name, res, truth = truth_mat, shared_cts = shared_celltypes, shared_samp = shared_samples, get_spill_cors = FALSE))) %>%
    dplyr::select(method, cors) %>%
    unnest(cols = c(cors)) %>%
    mutate(ref = refName,
           val = valDataset)

  return(out)


}, mc.cores = 10) %>%
  bind_rows()

all_cors_ref_combined <- all_cors %>%
  # filter(ref == "sc_pan_cancer") %>%
  #filter(method %in% c("xCell2", "dtangle")) %>%
  group_by(method, val) %>%
  dplyr::summarise(cors_list = list(cor),
                   n_ct_samples = list(n)) %>%
  rowwise() %>%
  mutate(n_ct_samples = list(ifelse(n_ct_samples > 30, 30, n_ct_samples))) %>%
  mutate(val_rho = list(combineRhos(rhos = cors_list, sample_sizes = log(n_ct_samples), use_median = FALSE, summarize = FALSE)),
         n_val_cts = length(cors_list))


data_combined <- all_cors_ref_combined %>%
  mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
  unnest(val_rho) %>%
  ungroup() %>%
  mutate(method = factor(method),
         val = factor(val)) %>%
  select(-c(cors_list, n_ct_samples))

train_ds <- data_combined %>%
  group_by(is_xcell2, method, val) %>%
  summarise(val_rho = mean(val_rho)) %>%
  group_by(is_xcell2, val) %>%
  top_n(1, wt=val_rho) %>%
  group_by(val) %>%
  mutate(delta_val_rho = val_rho[is_xcell2 == "yes"] - val_rho[is_xcell2 == "no"]) %>%
  select(val, delta_val_rho) %>%
  ungroup() %>%
  unique() %>%
  arrange(delta_val_rho) %>%
  pull(val)

# all
val       delta_val_rho
<fct>             <dbl>
1 GSE93722       -0.105
2 GSE115823      -0.103
3 GSE130824      -0.0603
4 GSE120444      -0.0519
5 GSE127813      -0.0379
6 BG_blood       -0.0169
7 GSE107011      -0.0154
8 GSE106898      -0.0152
9 GSE77343       -0.0149
10 GSE65133       -0.00364

# kass_blood
val       delta_val_rho
<fct>             <dbl>
1 GSE77343       -0.0586
2 GSE130824      -0.0491
3 GSE20300       -0.0130
4 GSE127813      -0.00838
5 GSE107572      -0.00499

# ts_blood
# A tibble: 16 × 2
val       delta_val_rho
<fct>             <dbl>
1 GSE77344       -0.239
2 GSE77343       -0.0987
3 GSE127813      -0.0820
4 BG_blood       -0.0652
5 GSE106898      -0.0364
6 GSE107011      -0.0257
7 GSE130824      -0.0173

# lm22
val               delta_val_rho
<fct>                     <dbl>
1 GSE115823              -0.163
2 BG_blood               -0.134
3 GSE107011              -0.119
4 GSE20300               -0.103
5 GSE106898              -0.0730
6 GSE107990              -0.0554
7 GSE127813              -0.0404
8 GSE130824              -0.0328
9 GSE65133               -0.0140
10 SDY420                 -0.00773

GSE64655

# train_ds <- c("GSE115823", "GSE20300", "GSE130824", "GSE65133", "GSE127813","GSE107011")  # new best

GSE115823 - Nasal_Asthma
GSE20300 - blood count
GSE64655 -C ell Type Purification
# train_ds <- c("GSE115823", "GSE20300", "GSE64655")  # new best


# Run tuning --------------

train_ds <- unique(ref_val_pairs$val_dataset)

vals.refs.res.sub <- vals.refs.res %>%
  filter(val_dataset %in% train_ds)

ref_val_pairs.sub <- ref_val_pairs %>%
  filter(val_dataset %in% train_ds)



xCell2results <- parallel::mclapply(1:nrow(vals.refs.res.sub), function(i){


      print(paste0("-------------------- ", i, "/", nrow(vals.refs.res.sub), " --------------------"))

      # Load data
      val_ref <- paste0(vals.refs.res.sub[i,]$val_dataset, "_", vals.refs.res.sub[i,]$ref_name[[1]])
      print(val_ref)
      mix.in <- cyto.vals$mixtures[[vals.refs.res.sub[i,]$val_type]][[vals.refs.res.sub[i,]$val_dataset]]
      ref.in <- refsRDSList[[vals.refs.res.sub[i,]$ref_type]][[vals.refs.res.sub[i,]$ref_name[[1]]]]$ref
      labels <- refsRDSList[[vals.refs.res.sub[i,]$ref_type]][[vals.refs.res.sub[i,]$ref_name[[1]]]]$labels
      lineage_file <- refsRDSList[[vals.refs.res.sub[i,]$ref_type]][[vals.refs.res.sub[i,]$ref_name[[1]]]]$lineage_file
      refType <- ifelse(vals.refs.res.sub[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res.sub[i,]$ref_type)
      valType <- vals.refs.res.sub[i,]$val_type
      valDataset <- vals.refs.res.sub[i,]$val_dataset
      refName <- vals.refs.res.sub[i,]$ref_name
      valDataType <- vals.refs.res.sub[i,]$val_data_type[[1]]

      file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/", val_ref, ".xcell2object.rds")

      if (file.exists(file)) {
        tryCatch({
          xcell2_object <- readRDS(file)
        }, error = function(e) {
          # If the specific error occurs, delete the file and print a message
          if (grepl("error reading from connection", e$message)) {
            file.remove(file)
            message(paste("File", file, "has been deleted due to error reading from connection."))
          } else {
            stop(e)  # Re-throw other errors
          }
        })
      }

      if (file.exists(file)) {
        print(paste0("Loading ", val_ref, " xCell2 Object.."))
        xcell2_object <- readRDS(file)
      }else{

        xcell2_object <- xCell2::xCell2Train(ref = ref.in, labels = labels, mix = mix.in, ref_type = refType, lineage_file = lineage_file,
                                             num_threads = 35, return_analysis = FALSE, return_signatures = TRUE,
                                             use_ontology = TRUE,
                                             probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.333, 0.4, 0.49),
                                             diff_vals = round(c(log(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20), log2(30), log2(40), log2(50), log2(100)), 3),
                                             min_frac_ct_passed = 0.25,
                                             min_genes = 3,
                                             max_genes = 500
        )
        saveRDS(xcell2_object, file)
      }

      signatures <- xcell2_object@signatures
      truth_mat <- cyto.vals$truth[[valType]][[valDataset]]
      sigs_celltypes <- unique(unlist(lapply(names(xcell2_object@signatures), function(x){strsplit(x, "#")[[1]][1]})))
      celltypes <- intersect(sigs_celltypes, rownames(truth_mat))

      res <- lapply(celltypes, function(ctoi){
        signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
        samples <- intersect(colnames(truth_mat), colnames(mix.in))

        mix_ranked <- singscore::rankGenes(mix.in[,samples])
        scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
          suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
        })

        fracs <- truth_mat[ctoi, samples]
        c <- apply(scores, 2, function(x){
          cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
        })

        enframe(c, name = "sig_name", value = "cor") %>%
          separate(sig_name, into = c("ct", "sig_prob", "sig_diff", "n_sig_genes", "frac_ct_passed"), sep = "_", convert = TRUE) %>%
          return(.)
      })

      res <- bind_rows(res)
      return(res)

}, mc.cores = 6)
saveRDS(xCell2results, "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/xcell2object.rds")

# Analyze tuning results ----------------------

xCell2results <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/xcell2object.rds")

vals.refs <- vals.refs.res.sub %>% select(val_dataset, ref_name, n_val_samples, ref_type, val_data_type)

vals.refs$res <- xCell2results



res.in <- vals.refs %>%
  unnest(res) %>%
  filter(sig_diff <= log2(101)) %>%
  mutate(cor = ifelse(is.na(cor), 0, cor)) %>%
  mutate(n_val_samples = ifelse(n_val_samples > 30, 30, n_val_samples)) %>%
  mutate(cor = ifelse(cor == 1, 0.99999, cor)) %>%
  mutate(cor = ifelse(cor == -1, -0.99999, cor)) %>%
  mutate(z = 0.5 * log((1 + cor) / (1 - cor))) %>%
  mutate(z_weighted = log(n_val_samples) * z) %>%
  mutate(rho_weighted = (exp(2 * z_weighted) - 1) / (exp(2 * z_weighted) + 1))


data <- res.in %>%
  filter(ref_name == "sc_pan_cancer")

data$val_dataset <- as.factor(data$val_dataset)
data$ref_name <- as.factor(data$ref_name)
data$celltype <- as.factor(data$ct)
data$ref_type <- as.factor(data$ref_type)
data$val_data_type <- as.factor(data$val_data_type)

# Create bins for each parameter
data <- data %>%
  mutate(
    n_sig_genes_bin = cut(n_sig_genes, breaks = c(3, 5, seq(10, max(n_sig_genes, na.rm = TRUE), by = 10)), include.lowest = TRUE),
    sig_prob = factor(sig_prob),
    sig_diff = factor(sig_diff),
    ref_type = factor(ref_type),
    val_data_type = factor(val_data_type),
    frac_ct_passed_bin = cut(frac_ct_passed, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)
  )



# Fit the model including the interaction between sig_prob and sig_diff
model_interaction_prob_diff_data <- lme4::lmer(rho_weighted ~ sig_prob * sig_diff + (1 | val_dataset) + (1 | ref_name) + (1 | ct), data = data)
# model_interaction_prob_diff_data <- lme4::lmer(rho_weighted ~ sig_prob * sig_diff + (1 | val_dataset) + (1 | ct), data = data)


# Summarize the model
summary(model_interaction_prob_diff_data)

# Get the predicted values for the interaction between sig_prob, sig_diff, ref_type, and val_data_type
effects_prob_diff_data <- ggeffects::ggpredict(model_interaction_prob_diff_data, terms = c("sig_prob", "sig_diff"))

# Convert the effects to a data frame for ggplot2
effects_df_prob_diff_data <- as.data.frame(effects_prob_diff_data)
colnames(effects_df_prob_diff_data)[1] <- "sig_prob"
colnames(effects_df_prob_diff_data)[6] <- "sig_diff"

# Define the custom color palette
color_palette <-c(colorRampPalette(c("darkred", "#FF6A6A", "#E0EEEE",  "royalblue1", "darkblue"))(50))


# Create a heatmap with custom color palette for each ref_type and val_data_type
ggplot(effects_df_prob_diff_data, aes(x = sig_prob, y = sig_diff, fill = predicted)) +
  geom_tile() +
  scale_fill_gradientn(colors = color_palette) +
  labs(title = "Interaction Effect of sig_prob and sig_diff Parameters",
       x = "sig_prob",
       y = "sig_diff",
       fill = "Predicted\nRho Weighted") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


model_interaction_genes_frac <- lme4::lmer(rho_weighted ~ n_sig_genes_bin * frac_ct_passed_bin + (1 | val_dataset) + (1 | ref_name) + (1 | celltype), data = data)
# model_interaction_genes_frac <- lme4::lmer(rho_weighted ~ n_sig_genes_bin * frac_ct_passed_bin + (1 | val_dataset) + (1 | celltype), data = data)


# Summarize the model
summary(model_interaction_genes_frac)

# Get the predicted values for the interaction between n_sig_genes_bin and frac_ct_passed_bin
effects_genes_frac <- ggeffects::ggpredict(model_interaction_genes_frac, terms = c("n_sig_genes_bin", "frac_ct_passed_bin"))

# Convert the effects to a data frame for ggplot2
effects_df_genes_frac <- as.data.frame(effects_genes_frac)
colnames(effects_df_genes_frac)[1] <- "n_sig_genes_bin"
colnames(effects_df_genes_frac)[6] <- "frac_ct_passed_bin"



# Define the custom color palette
#color_palette <- colorRampPalette(c("darkred", "red", "white", "lightblue1", "darkblue"))(50)

# Create a heatmap with custom color palette
ggplot(effects_df_genes_frac, aes(x = n_sig_genes_bin, y = frac_ct_passed_bin, fill = predicted)) +
  geom_tile() +
  scale_fill_gradientn(colors = color_palette) +
  labs(title = "Interaction Effect of n_sig_genes_bin and frac_ct_passed_bin on rho_weighted",
       x = "n_sig_genes_bin",
       y = "frac_ct_passed_bin",
       fill = "Predicted rho_weighted") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

