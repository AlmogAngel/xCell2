library(tidyverse)
library(parallel)


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
      select(celltype, cor, p_value, n) %>%
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
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

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
  select(-n) %>%
  ungroup()


# Load results from run_sig_gen_tuning.R
xcell2_resfiles <- list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/sig_gen_tuning_3", full.names = TRUE, pattern = "*.rds")
for (xcell2_resfile in xcell2_resfiles) {

  y <- readRDS(xcell2_resfile)
  y <- lapply(y, function(x){round(x, 3)})
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
    group_by(method, ref) %>%
    dplyr::summarise(cors_list = list(cor),
                     n_ct_samples = list(n)) %>% # Rhos are weighted by number of samples per cell type
    rowwise() %>%
    mutate(ref_rho = list(combineRhos(rhos = cors_list, sample_sizes = log(n_ct_samples), use_median = FALSE, summarize = FALSE)),
           n_val_cts = length(cors_list))

  data_combined <- all_cors_ref_combined %>%
    mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
    unnest(ref_rho) %>%
    ungroup() %>%
    mutate(method = factor(method),
           ref = factor(ref)) %>%
    select(-c(cors_list, n_ct_samples))

  p <- ggplot(data_combined, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
    geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
    scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
    theme_minimal() +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) +
    coord_flip() +
    facet_wrap(~ ref, scales = "free", ncol = 1) +
    tidytext::scale_x_reordered() +
    scale_x_discrete(labels=bold_xCell2_labels) +
    labs(x="", y="Weighted Spearman Rho", title = gsub("_res.rds", "", basename(xcell2_resfile))) +
    guides(fill=FALSE)

  ggsave(filename = paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/sig_gen_tuning_3/plots/", gsub("_res.rds", "", basename(xcell2_resfile)), ".pdf"), plot = p, width = 441, height = 416, dpi = 300, units = "mm", device = "pdf")


}






