library(tidyverse)


round_results <- 3

# Spillover analysis

getCors <- function(ref, res, truth, shared_cts, shared_samp, get_spill_cors, cor_method){

  makeGEPMat <- function(ref, labels){

    celltypes <- unique(labels$label)

    gep_mat <- sapply(celltypes, function(type){
      type_samples <- labels[,2] == type
      if (sum(type_samples) == 1) {
        type_vec <- as.vector(ref[,type_samples])
      }else{
        type_vec <- Rfast::rowMedians(as.matrix(ref[,type_samples]))
      }
    })
    rownames(gep_mat) <- rownames(ref)

    return(gep_mat)
  }
  getCellTypeCorrelation <- function(gep_mat, ref_type){

    celltypes <- colnames(gep_mat)

    if (ref_type != "sc") {

      # Use top 10% most variable genes
      genes_var <- apply(gep_mat, 1, var)
      most_var_genes_cutoff <- quantile(genes_var, 0.9, na.rm=TRUE)
      gep_mat <- gep_mat[genes_var > most_var_genes_cutoff,]

    }else{

      # Use top 1% most variable genes
      genes_var <- apply(gep_mat, 1, var)
      most_var_genes_cutoff <- quantile(genes_var, 0.99, na.rm=TRUE)
      gep_mat <- gep_mat[genes_var > most_var_genes_cutoff,]

    }


    # Make correlation matrix
    cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
    lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

    # TODO: Change for loop to apply function to measure time
    for (i in 1:nrow(lower_tri_coord)) {
      celltype_i <- rownames(cor_mat)[lower_tri_coord[i, 1]]
      celltype_j <- colnames(cor_mat)[lower_tri_coord[i, 2]]
      cor_mat[lower_tri_coord[i, 1], lower_tri_coord[i, 2]] <- cor(gep_mat[,celltype_i], gep_mat[,celltype_j], method = "spearman")
      cor_mat[lower_tri_coord[i, 2], lower_tri_coord[i, 1]] <- cor(gep_mat[,celltype_i], gep_mat[,celltype_j], method = "spearman")
    }

    return(cor_mat)
  }


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

    #     source("/bigdata/almogangel/xCell2/dev_scripts/xCell2Train_tmp_version.R")
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
    dep_list <- getDependencies(ref.in$lineage_file)


    # # Generate all pairs of cell types
    # celltype_pairs <- crossing(predicted_celltype = unique(df$celltype), true_celltype = unique(df$celltype))
    # celltype_pairs <- celltype_pairs %>%
    #   rowwise() %>%
    #   mutate(true_celltype = ifelse(true_celltype %in% c(predicted_celltype, unlist(dep_list[predicted_celltype])), NA, true_celltype)) %>%
    #   drop_na()


    # Generate pairs of CTOI vs. most similar cell type
    gep_mat <- makeGEPMat(ref.in$ref, ref.in$labels)
    ref_type <- ifelse(ref %in% c("ts_blood", "sc_pan_cancer"), "sc", "rnaseq")
    cor_mat <- getCellTypeCorrelation(gep_mat, ref_type)
    most_similar_cts <- sapply(celltypes, function(ct){
      celltypes2use <-  celltypes[!celltypes %in% c(ct, unlist(dep_list[ct]))]
      names(sort(cor_mat[ct, celltypes2use], decreasing = TRUE))[1]
    })
    celltype_pairs <- tibble(predicted_celltype = most_similar_cts, true_celltype = celltypes)


    cors.out <- celltype_pairs %>%
      rowwise() %>%
      mutate(cor = cor(
        df %>% filter(celltype == predicted_celltype) %>% pull(prediction),
        df %>% filter(celltype == true_celltype) %>% pull(truth),
        method = cor_method, use = "pairwise.complete.obs"),
        p_value = cor.test(
          df %>% filter(celltype == predicted_celltype) %>% pull(prediction),
          df %>% filter(celltype == true_celltype) %>% pull(truth),
          method = cor_method, exact = FALSE)$p.value,
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
        cor = cor(truth, prediction, method = cor_method, use = "pairwise.complete.obs"),
        p_value = cor.test(truth, prediction, method = cor_method, exact = FALSE)$p.value,
        n = sum(!is.na(truth) & !is.na(prediction))
      ) %>%
      mutate(cor = ifelse(is.na(cor), 0, cor)) %>%
      return(.)
  }

}
getDependencies <- function(lineage_file_checked){
  ont <- read_tsv(lineage_file_checked, show_col_types = FALSE) %>%
    mutate_all(as.character)

  celltypes <- pull(ont[,2])
  celltypes <- gsub("_", "-", celltypes)
  dep_list <- vector(mode = "list", length = length(celltypes))
  names(dep_list) <- celltypes

  for (i in 1:nrow(ont)) {
    descendants <-  gsub("_", "-", strsplit(pull(ont[i,3]), ";")[[1]])
    descendants <- descendants[!is.na(descendants)]

    ancestors <-  gsub("_", "-", strsplit(pull(ont[i,4]), ";")[[1]])
    ancestors <- ancestors[!is.na(ancestors)]

    dep_list[[i]] <- list("descendants" = descendants, "ancestors" = ancestors)

  }

  return(dep_list)
}
combineRhos <- function(rhos, sample_sizes = NULL, use_median = FALSE, summarize = FALSE){

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
makeGEPMat <- function(ref, labels){

  celltypes <- unique(labels$label)

  gep_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      type_vec <- Rfast::rowMedians(as.matrix(ref[,type_samples]))
    }
  })
  rownames(gep_mat) <- rownames(ref)

  return(gep_mat)
}
getCellTypeCorrelation <- function(gep_mat, ref_type){

  celltypes <- colnames(gep_mat)

  if (ref_type != "sc") {

    # Use top 10% most variable genes
    genes_var <- apply(gep_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.9, na.rm=TRUE)
    gep_mat <- gep_mat[genes_var > most_var_genes_cutoff,]

  }else{

    # Use top 1% most variable genes
    genes_var <- apply(gep_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.99, na.rm=TRUE)
    gep_mat <- gep_mat[genes_var > most_var_genes_cutoff,]

  }


  # Make correlation matrix
  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

  # TODO: Change for loop to apply function to measure time
  for (i in 1:nrow(lower_tri_coord)) {
    celltype_i <- rownames(cor_mat)[lower_tri_coord[i, 1]]
    celltype_j <- colnames(cor_mat)[lower_tri_coord[i, 2]]
    cor_mat[lower_tri_coord[i, 1], lower_tri_coord[i, 2]] <- cor(gep_mat[,celltype_i], gep_mat[,celltype_j], method = "spearman")
    cor_mat[lower_tri_coord[i, 2], lower_tri_coord[i, 1]] <- cor(gep_mat[,celltype_i], gep_mat[,celltype_j], method = "spearman")
  }

  return(cor_mat)
}


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


files <- list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/", pattern = "*spillAlpha*", full.names = TRUE)
files <- c("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params//xcell2_benchmark_results_spillAlpha_0.rds",
           "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params//xcell2_benchmark_results_spillAlpha_0.25.rds",
           "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params//xcell2_benchmark_results_spillAlpha_0.5.rds",
           "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params//xcell2_benchmark_results_spillAlpha_0.75.rds",
           "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params//xcell2_benchmark_results_spillAlpha_1.rds")
use_partial_cors <- FALSE

spill_cor_res_table <- lapply(files, function(file){

  print(file)
  vals.refs.res$res <- readRDS(file)

  parallel::mclapply(1:nrow(vals.refs.res), function(i){


    valType <- vals.refs.res[i,]$val_type
    valDataset <- vals.refs.res[i,]$val_dataset
    refName <- vals.refs.res[i,]$ref_name
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))
    dep_list <- getDependencies(ref.in$lineage_file)
    gep_mat <- makeGEPMat(ref.in$ref, ref.in$labels)
    ref_type <- ifelse(refName %in% c("ts_blood", "sc_pan_cancer"), "sc", "rnaseq")
    cor_mat <- getCellTypeCorrelation(gep_mat, ref_type)
    truth <- cyto.vals$truth[[valType]][[valDataset]]
    truth <- round(truth, round_results)
    res <- vals.refs.res[i,]$res[[1]]
    res <- round(res, round_results)

    celltypes <- intersect(rownames(res), rownames(truth))
    samples <- intersect(colnames(res), colnames(truth))

    if (use_partial_cors) {


      spill_cor_res <- lapply(celltypes, function(ct){

       celltypes2use <-  celltypes[!celltypes %in% c(unlist(dep_list[ct]))]

        truth_tmp <- t(truth[celltypes2use,])
        y_truth <- truth_tmp[,ct]
        y_truth <- y_truth[!is.na(y_truth)]
        x_scores <- res[ct, ]

        samples2use <- intersect(names(x_scores), names(y_truth))
        x_scores <- x_scores[samples2use]
        y_truth <- y_truth[samples2use]


        z_table <- truth_tmp[names(y_truth),]
        z_table <- z_table[,colnames(z_table) != ct]
        z_table <- z_table[,colSums(is.na(z_table)) == 0]

        cor_table <- cbind("X" = x_scores, "Y" = y_truth, z_table)
        partial_corr_complete <- ppcor::pcor(cor_table, method = "pearson")
        pcor_res <- partial_corr_complete$estimate[,1]
        names(pcor_res) <- colnames(cor_table)

        true_cor <- cor(x_scores, y_truth, method = "pearson")

        tibble(sig_ct = ct, direct_cor = true_cor, partial_cor = pcor_res["Y"], n_samples = length(x_scores))

      }) %>%
        bind_rows()

      y <- vals.refs.res[i,] %>%
        select(method:shared_celltypes)

      spill_cor_res <- cbind(y, spill_cor_res) %>%
        return(.)


    }else{
      ct2most_simillar <- sapply(celltypes, function(ct){
        celltypes2use <-  celltypes[!celltypes %in% c(ct, unlist(dep_list[ct]))]
        names(sort(cor_mat[ct, celltypes2use], decreasing = TRUE))[1]
      })

      y <- vals.refs.res[i,] %>%
        select(method:shared_celltypes)

      spill_cor_res <- tibble(sig_ct = celltypes, most_sim_truth_ct = ct2most_simillar)
      spill_cor_res <- cbind(y, spill_cor_res)

      spill_cor_res %>%
        rowwise() %>%
        mutate(spill_cor = cor(res[sig_ct, samples], truth[most_sim_truth_ct, samples], method = "pearson", use = "pairwise.complete.obs")) %>%
        mutate(true_cor = cor(res[sig_ct, samples], truth[sig_ct, samples], method = "pearson", use = "pairwise.complete.obs")) %>%
        ungroup() %>%
        return(.)
    }


  }, mc.cores = 20) %>%
    bind_rows() %>%
    mutate(spill_alpha = gsub(".rds", "", gsub("xcell2_benchmark_results_spillAlpha_", "", basename(file),)))

}) %>%
  bind_rows()

# spill_cor_res_table_adjusted <- as_tibble(spill_cor_res_table) %>%
#   group_by(method, spill_alpha, ref_name, val_dataset) %>%
#   dplyr::summarise(direct_cor_list = list(direct_cor),
#                    partial_cor_list = list(partial_cor),
#                    n_ct_samples = list(n_samples)) %>% # Rhos are weighted by number of samples per cell type
#   rowwise() %>%
#   mutate(n_ct_samples = list(ifelse(n_ct_samples > 30, 30, n_ct_samples))) %>%
#   mutate(direct_cor_adjusted = list(combineRhos(rhos = direct_cor_list, sample_sizes = log(n_ct_samples), use_median = FALSE, summarize = FALSE)),
#          partial_cor_adjusted = list(combineRhos(rhos = partial_cor_list, sample_sizes = log(n_ct_samples), use_median = FALSE, summarize = FALSE)),
#          n_val_cts = length(direct_cor_list))
#
#
# d_table <- spill_cor_res_table_adjusted %>%
#   unnest_longer(direct_cor_adjusted) %>%
#   select(method, spill_alpha, ref_name, val_dataset, direct_cor_adjusted)
#
# p_table <- spill_cor_res_table_adjusted %>%
#   unnest_longer(partial_cor_adjusted) %>%
#   select(method, spill_alpha, ref_name, val_dataset, partial_cor_adjusted)
#
#
# final_partial_adjusted_cors <- cbind(d_table, "partial_cor_adjusted" = p_table$partial_cor_adjusted)
#
#
# as_tibble(final_partial_adjusted_cors) %>%
#   drop_na() %>%
#   pivot_longer(cols = c(direct_cor_adjusted, partial_cor_adjusted), names_to = "cor_type", values_to = "cor") %>%
#   ggplot(., aes(x=spill_alpha, y=cor, fill=cor_type)) +
#   #coord_flip(ylim = c(0, 0.8)) +
#   geom_boxplot()



spill_cor_res_table <- spill_cor_res_table %>%
  pivot_longer(cols = c(spill_cor, true_cor), names_to = "cor_type", values_to = "cor")

spill_cor_res_table[is.na(spill_cor_res_table$cor),]$cor <- 0

spill_cor_res_table$cor_type <- factor(spill_cor_res_table$cor_type, levels = c("true_cor", "spill_cor"))
spill_cor_res_table %>%
  #filter(ref_type != "sc") %>%
  ggplot(., aes(x=spill_alpha, y=cor, fill=cor_type)) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("spill_cor"="tomato", "true_cor"="#7AC5CD"),
                    labels = c("Direct", "Spill")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "bottom") +
  labs(x="Alpha", y="Pearson r", fill="Correlation type")




spill_cor_res_table %>%
  filter(ref_name == "sc_pan_cancer"  & cor_type == "true_cor") %>%
  filter(spill_alpha %in% c(0, 0.5)) %>%
  pivot_wider(names_from = spill_alpha, values_from = cor) %>%
  select(ref_name, val_dataset, ref_type, sig_ct, most_sim_truth_ct, cor_type, `0`, `0.5`) %>%
  mutate(delta = `0` - `0.5`) %>%
  arrange(delta) %>%
  print(n=Inf)



