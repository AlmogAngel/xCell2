library(tidyverse)

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




# Choose 8 datasets that xCell 2.0 preforms worst (remove) -----------------

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



y <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/sig_filt_tuning_2/grid_36_res.rds")
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
  group_by(method, val) %>%
  dplyr::summarise(cors_list = list(cor),
                   n_ct_samples = list(n)) %>%
  rowwise() %>%
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



# 1 GSE64655      -0.397   - Bulk RNA-seq	Cell Type Purification	Blood
# 2 GSE93722      -0.305   - Bulk RNA-seq	Flow Cytometry	Lymph
# 3 GSE20300      -0.144   - Bulk MicroArray	Blood Count	Blood
# 4 GSE77344      -0.0674  - Bulk MicroArray	Flow Cytometry	Blood
# 5 GSE77343      -0.0635  - Bulk MicroArray	Flow Cytometry	Blood
# 6 GSE106898      0.000552  - Bulk MicroArray	Flow Cytometry	Blood

train_ds <- c("GSE64655", "GSE93722", "GSE20300", "GSE77344", "GSE77343", "GSE106898")


# Make parameters grid -----------------

# Best parameters (all validations):
# min_genes = 3
# max_genes = 200
# min_frac_ct_passed = 0.5
# probs = c(0.010, 0.050, 0.100, 0.150, 0.200, 0.333, 0.400)
# diff_vals = round(c(log(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3)

# Tuning for signature generation:
# use_ontology = c(TRUE, FALSE)
# min_genes = c(3, 5, 10,)
# max_genes = c(100, 200, 300)
# probs = list(c(0.01, 0.05, 0.1, 0.15, 0.2, 0.333, 0.4),
#              c(0.05, 0.1, 0.15, 0.2, 0.333, 0.4),
#              c(0.01, 0.05, 0.1, 0.15, 0.2, 0.333))
# diff_vals = list(round(c(log(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3),
#                  round(c(log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5),log2(10), log2(20)), 3),
#                  round(c(log2(2), log2(2.5), log2(3), log2(4), log2(5),log2(10), log2(20)), 3))
# min_frac_ct_passed = c(0.25, 0.5, 0.75, 0.9)
# grid <- expand_grid(use_ontology, min_genes, max_genes, probs, diff_vals, min_frac_ct_passed)

# Best parameters for signature generation:
# grid[which.max(grid_res),]
# use_ontology_best = TRUE
# min_genes_best = 10
# max_genes_best = 300
# probs_best = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.333)
# diff_vals_best = round(c(log(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3)
# min_frac_ct_passed_best = 0.75
#
#
# # Tuning for spillover:
# top_spill_value = c(0.1, 0.25, 0.5)
# sc_spill_relaxing_factor = c(1, 0.5, 0.25)
# spillover_alpha = c(0.25, 0.5, 0.75)
# grid <- expand_grid(top_spill_value, sc_spill_relaxing_factor, spillover_alpha)

# Run tuning --------------

vals.refs.res.sub <- vals.refs.res %>%
  filter(val_dataset %in% train_ds)

ref_val_pairs.sub <- ref_val_pairs %>%
  filter(val_dataset %in% train_ds)


grid_res <- c()
for (j in 1:nrow(grid)) {

  file <-  paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_train_tuning/spillover/grid_", j, "_res.rds")


  if (file.exists(file)) {
    print(paste0("Gird: ", j, " exist..."))
    xCell2results <- readRDS(file)
  }else{
    print(paste0("Gird: ", j, "/", nrow(grid)))
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

      file2 <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_train_tuning/spillover/grid_", j, "##", val_ref, ".xcell2object.rds")


      if (file.exists(file2)) {
        print(paste0("Loading ", val_ref, " xCell2 Object.."))
        xcell2_object <- readRDS(file2)
      }else{
        # xcell2_object <- xCell2::xCell2Train(ref = ref.in, labels = labels, mix = mix.in, filtering_data = NULL, ref_type = refType, lineage_file = lineage_file,
        #                                      num_threads = 30, return_analysis = FALSE, return_signatures = TRUE, filter_sigs = FALSE,
        #                                      use_ontology = grid$use_ontology[j],
        #                                      probs = grid$probs[j][[1]],
        #                                      diff_vals = grid$diff_vals[j][[1]],
        #                                      min_frac_ct_passed = grid$min_frac_ct_passed[j],
        #                                      min_genes = grid$min_genes[j],
        #                                      max_genes = grid$max_genes[j]
        # )
        xcell2_object <- xCell2::xCell2Train(ref = ref.in, labels = labels, mix = mix.in, filtering_data = NULL, ref_type = refType, lineage_file = lineage_file,
                                             num_threads = 35, return_analysis = FALSE, return_signatures = FALSE, filter_sigs = FALSE,
                                             use_ontology = use_ontology_best,
                                             probs = probs_best,
                                             diff_vals = diff_vals_best,
                                             min_frac_ct_passed = min_frac_ct_passed_best,
                                             min_genes = min_genes_best,
                                             max_genes = max_genes_best,
                                             top_spill_value = grid[j,]$top_spill_value,
                                             sc_spill_relaxing_factor = grid[j,]$sc_spill_relaxing_factor
        )
        saveRDS(xcell2_object, file2)
      }


      # res <- xCell2::xCell2Analysis(mix.in, xcell2object = xcell2_object, raw_scores = TRUE,
      #                               spillover = FALSE, spillover_alpha = 0, num_threads = 20
      #                               )
      res <- xCell2::xCell2Analysis(mix.in, xcell2object = xcell2_object, raw_scores = FALSE,
                                    spillover = TRUE, spillover_alpha = grid[j,]$spillover_alpha, num_threads = 20
      )

      return(res)

    }, mc.cores = 6)
    saveRDS(xCell2results, file)
  }


  xCell2results <- lapply(xCell2results, function(x){round(x, 3)})
  vals.refs.res.sub$res <- xCell2results

  all_cors <- parallel::mclapply(1:nrow(ref_val_pairs.sub), function(i){


    valType <- ref_val_pairs.sub[i,]$val_type
    valDataset <- ref_val_pairs.sub[i,]$val_dataset
    refName <- ref_val_pairs.sub[i,]$ref_name
    truth_mat <- cyto.vals$truth[[valType]][[valDataset]]
    truth_mat <- round(truth_mat, 3)

    cyto.Res.tmp <- cyto.Res %>%
      filter(ref_name == refName & val_dataset == valDataset)

    yy <- cyto.Res.tmp[cyto.Res.tmp$method == "BayesPrism",]
    yy$method <- "xCell2"
    yy$res <- pull(filter(vals.refs.res.sub, ref_name == refName & val_dataset == valDataset), res)
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


  }, mc.cores = 15) %>%
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
    dplyr::select(-c(cors_list, n_ct_samples))


  mean_res <- data_combined %>%
    filter(method == "xCell2") %>%
    group_by(ref) %>%
    summarise(m_rho = median(ref_rho)) %>%
    pull(m_rho) %>%
    mean()


  grid_res <- c(grid_res, mean_res)

}

saveRDS(grid_res, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_train_tuning/spillover/res_1.rds")


# Get results with the tuned parameters ----------



if (FALSE) {
  best_grid_id <- which.max(grid_res)
  best_grid <- grid[best_grid_id,]


  vals.refs.res.sub2 <- vals.refs.res %>%
    filter(!val_dataset %in% train_ds)

  ref_val_pairs.sub2 <- ref_val_pairs %>%
    filter(!val_dataset %in% train_ds)


  xCell2results <- parallel::mclapply(1:nrow(vals.refs.res.sub2), function(i){


    print(paste0("-------------------- ", i, "/", nrow(vals.refs.res.sub2), " --------------------"))

    # Load data
    val_ref <- paste0(vals.refs.res.sub2[i,]$val_dataset, "_", vals.refs.res.sub2[i,]$ref_name[[1]])
    print(val_ref)
    mix.in <- cyto.vals$mixtures[[vals.refs.res.sub2[i,]$val_type]][[vals.refs.res.sub2[i,]$val_dataset]]
    ref.in <- refsRDSList[[vals.refs.res.sub2[i,]$ref_type]][[vals.refs.res.sub2[i,]$ref_name[[1]]]]$ref
    labels <- refsRDSList[[vals.refs.res.sub2[i,]$ref_type]][[vals.refs.res.sub2[i,]$ref_name[[1]]]]$labels
    lineage_file <- refsRDSList[[vals.refs.res.sub2[i,]$ref_type]][[vals.refs.res.sub2[i,]$ref_name[[1]]]]$lineage_file
    refType <- ifelse(vals.refs.res.sub2[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res.sub2[i,]$ref_type)
    valType <- vals.refs.res.sub2[i,]$val_type
    valDataset <- vals.refs.res.sub2[i,]$val_dataset
    refName <- vals.refs.res.sub2[i,]$ref_name
    valDataType <- vals.refs.res.sub2[i,]$val_data_type[[1]]

    file_best <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_train_tuning/best_grid_", best_grid_id, "##", val_ref, ".xcell2object.rds")


    if (file.exists(file_best)) {
      print(paste0("Loading ", val_ref, " xCell2 Object.."))
      xcell2_object <- readRDS(file_best)
    }else{

      xcell2_object <- xCell2::xCell2Train(ref = ref.in, labels = labels, mix = mix.in, filtering_data = NULL, ref_type = refType, lineage_file = lineage_file,
                                           num_threads = 25, return_analysis = FALSE, return_signatures = FALSE, filter_sigs = FALSE,
                                           use_ontology = use_ontology_best,
                                           probs = probs_best,
                                           diff_vals = diff_vals_best,
                                           min_frac_ct_passed = min_frac_ct_passed_best,
                                           min_genes = min_genes_best,
                                           max_genes = max_genes_best,
                                           top_spill_value = 0.5,
                                           sc_spill_relaxing_factor = 1
      )

      saveRDS(xcell2_object, file_best)
    }


    res <- xCell2::xCell2Analysis(mix.in, xcell2object = xcell2_object, raw_scores = FALSE,
                                  spillover = TRUE, spillover_alpha = 0.75, num_threads = 20
    )

    return(res)

  }, mc.cores = 8)
  saveRDS(xCell2results, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_train_tuning/best_grid_final.rds")

  y <- xCell2results
  y <- lapply(y, function(x){round(x, 3)})
  ref_val_pairs.sub2$res <- y

  all_cors <- parallel::mclapply(1:nrow(vals.refs.res.sub2), function(i){


    valType <- vals.refs.res.sub2[i,]$val_type
    valDataset <- vals.refs.res.sub2[i,]$val_dataset
    refName <- vals.refs.res.sub2[i,]$ref_name
    truth_mat <- cyto.vals$truth[[valType]][[valDataset]]


    cyto.Res.tmp <- cyto.Res %>%
      filter(ref_name == refName & val_dataset == valDataset)

    yy <- cyto.Res.tmp[cyto.Res.tmp$method == "BayesPrism",]
    yy$method <- "xCell2"
    yy$res <- pull(filter(vals.refs.res.sub2, ref_name == refName & val_dataset == valDataset), res)
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

  ggplot(data_combined, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
    geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
    scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
    theme_minimal() +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) +
    coord_flip() +
    facet_wrap(~ ref, scales = "free", ncol = 1) +
    tidytext::scale_x_reordered() +
    scale_x_discrete(labels=bold_xCell2_labels) +
    labs(x="", y="Weighted Spearman Rho") +
    guides(fill=FALSE)
}
