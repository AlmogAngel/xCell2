library(tidyverse)

use_mouse <- FALSE


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

if (use_mouse) {
  methods2use <- c("xCell2", "BayesPrism", "CIBERSORTx", "EPIC", "MCPcounter", "dtangle", "Bisque", "DWLS", "SCDC", "DeconRNASeq", "Scaden")
  refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val_mouse.rds")
  cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.mouse.rds")

  refList <- list(rna_seq = c(mixed = "igd", mixed = "mouse_rnaseq_data"),
                  array = c(),
                  sc = c(mixed = "tm_blood"))

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

  files <- list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/mouse/", pattern = ".cyto.res.rds", full.names = TRUE)
  cyto.Res <- lapply(files, function(f){

    f.in <- readRDS(f)

    if (colnames(f.in)[1] != "method") {
      m <- gsub(".cyto.res.rds", "", basename(f))
      f.in <- as_tibble(cbind(method = m, f.in))
    }


    f.in %>%
      dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
      group_by(across(-ncol(.))) %>%
      summarise(res = list(do.call(rbind, res)), .groups = 'drop')
  }) %>%
    do.call(rbind, .)

  cyto.Res[cyto.Res$method == "dwlr", ]$method <- "DWLR"
  cyto.Res[cyto.Res$method == "scdc", ]$method <- "SCDC"
  cyto.Res[cyto.Res$method == "bisque", ]$method <- "Bisque"

}else{
  methods2use <- c("xCell2", "BayesPrism", "CIBERSORTx", "DeconRNASeq", "EPIC", "MCPcounter", "dtangle", "Bisque", "DWLS", "SCDC", "Scaden")
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

    f.in <- readRDS(f)

    if (colnames(f.in)[1] != "method") {
      m <- gsub(".cyto.res.rds", "", basename(f))
      f.in <- as_tibble(cbind(method = m, f.in))
    }


    f.in %>%
      dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
      group_by(across(-ncol(.))) %>%
      summarise(res = list(do.call(rbind, res)), .groups = 'drop')
  }) %>%
    do.call(rbind, .)

  cyto.Res[cyto.Res$method == "dwlr", ]$method <- "DWLR"
  cyto.Res[cyto.Res$method == "scdc", ]$method <- "SCDC"
  cyto.Res[cyto.Res$method == "bisque", ]$method <- "Bisque"
}




get_xcell2_benchmark_results <- function(benchmark_table = vals.refs.res, vals2remove = c(), save_object = TRUE, dir = NULL, params, ncores, output_name){




  benchmark_table <- benchmark_table %>%
    filter(!val_dataset %in% vals2remove)

  print(paste0("Running benchmarking for ", length(unique(benchmark_table$val_dataset)), " validation datasets..."))

  xCell2results <- parallel::mclapply(1:nrow(benchmark_table), function(i){


    print(paste0("-------------------- ", i, "/", nrow(benchmark_table), " --------------------"))

    # Load data
    val_ref <- paste0(benchmark_table[i,]$val_dataset, "_", benchmark_table[i,]$ref_name[[1]])
    print(val_ref)
    mix.in <- cyto.vals$mixtures[[benchmark_table[i,]$val_type]][[benchmark_table[i,]$val_dataset]]
    ref.in <- refsRDSList[[benchmark_table[i,]$ref_type]][[benchmark_table[i,]$ref_name[[1]]]]$ref
    labels <- refsRDSList[[benchmark_table[i,]$ref_type]][[benchmark_table[i,]$ref_name[[1]]]]$labels
    lineage_file <- refsRDSList[[benchmark_table[i,]$ref_type]][[benchmark_table[i,]$ref_name[[1]]]]$lineage_file
    refType <- ifelse(benchmark_table[i,]$ref_type == "rna_seq", "rnaseq", benchmark_table[i,]$ref_type)
    valType <- benchmark_table[i,]$val_type
    valDataset <- benchmark_table[i,]$val_dataset
    refName <- benchmark_table[i,]$ref_name
    valDataType <- benchmark_table[i,]$val_data_type[[1]]


    if (save_object) {
      if (is.null(dir)) {
        errorCondition("Please provide dir parameter.")
      }
      file <- paste0(dir, "/", val_ref, ".xcell2object.rds")
    }

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
                                           human2mouse = params$human2mouse,
                                           num_threads = params$num_threads,
                                           min_pb_cells = params$min_pb_cells,
                                           min_pb_samples = params$min_pb_samples,
                                           min_sc_genes = params$min_sc_genes,
                                           use_ontology = params$use_ontology,
                                           return_signatures = params$return_signatures,
                                           return_analysis = params$return_analysis,
                                           use_sillover = params$use_sillover,
                                           spillover_alpha = params$spillover_alpha,
                                           top_spill_value = params$top_spill_value
                                           )

      saveRDS(xcell2_object, file)
    }


    if (!params$return_analysis) {
      # refIsSC <- ifelse(refType == "sc", TRUE, FALSE)
      refIsSC = FALSE
      res <- xCell2::xCell2Analysis(mix.in, xcell2object = xcell2_object, raw_scores = !params$use_sillover, ref_is_sc = refIsSC,
                                    spillover = params$use_sillover, spillover_alpha = params$spillover_alpha, num_threads = params$num_threads)
    }



    return(res)

  }, mc.cores = ncores)
  saveRDS(xCell2results, paste0(dir, "/", output_name))

  print(paste0("Done - xCell2 benchmarking analysis results located in: ", dir, "/", output_name))

}


get_xcell2_correlations <- function(benchmark_table = vals.refs.res, vals2remove = c(), xCell2results, round_results = 3, weight_rhos = TRUE,
                                    by_val = FALSE, ref2use = NULL, spillcors = FALSE, cMethod = "spearman"){

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


  cyto.Res <- cyto.Res %>%
    rowwise() %>%
    mutate(res = list(round(res, round_results)))


  if (spillcors) {

    spill_res <- parallel::mclapply(1:nrow(cyto.Res), function(i){


      valType <- cyto.Res[i,]$val_type
      valDataset <- cyto.Res[i,]$val_dataset
      refName <- cyto.Res[i,]$ref_name
      ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))
      dep_list <- getDependencies(ref.in$lineage_file)
      gep_mat <- makeGEPMat(ref.in$ref, ref.in$labels)
      ref_type <- ifelse(refName %in% c("ts_blood", "sc_pan_cancer"), "sc", "rnaseq")
      cor_mat <- getCellTypeCorrelation(gep_mat, ref_type)
      truth <- cyto.vals$truth[[valType]][[valDataset]]
      truth <- round(truth, roundResults)
      res <- cyto.Res[i,]$res[[1]]

      celltypes <- intersect(rownames(res), rownames(truth))
      samples <- intersect(colnames(res), colnames(truth))

      ct2most_simillar <- sapply(celltypes, function(ct){
        celltypes2use <-  celltypes[!celltypes %in% c(ct, unlist(dep_list[ct]))]
        names(sort(cor_mat[ct, celltypes2use], decreasing = TRUE))[1]
      })

      y <- cyto.Res[i,] %>%
        select(method:n_shared_celltypes)

      spill_cor_res <- tibble(sig_ct = celltypes, most_sim_truth_ct = ct2most_simillar)
      spill_cor_res <- cbind(y, spill_cor_res)

      spill_cor_res %>%
        rowwise() %>%
        mutate(spill_cor = cor(res[sig_ct, samples], truth[most_sim_truth_ct, samples], method = "pearson", use = "pairwise.complete.obs")) %>%
        mutate(direct_cor = cor(res[sig_ct, samples], truth[sig_ct, samples], method = "pearson", use = "pairwise.complete.obs")) %>%
        ungroup() %>%
        return(.)


    }, mc.cores = 20) %>%
      bind_rows()

      return(spill_res)

  }

  benchmark_table <- benchmark_table %>%
    filter(!val_dataset %in% vals2remove)

  if (!is.null(round_results)) {
    xCell2results <- lapply(xCell2results, function(x){round(x, round_results)})
  }

  cyto.Res <- cyto.Res %>%
    filter(!val_dataset %in% vals2remove)


  ref_val_pairs <- cyto.Res %>%
    group_by(ref_tissue, ref_name, val_type, val_dataset) %>%
    summarise(n = n()) %>%
    filter(n == length(methods2use)-1) %>%
    dplyr::select(-n) %>%
    ungroup()


  benchmark_table$res <- xCell2results


  all_cors <- parallel::mclapply(1:nrow(ref_val_pairs), function(i){


    valType <- ref_val_pairs[i,]$val_type
    valDataset <- ref_val_pairs[i,]$val_dataset
    refName <- ref_val_pairs[i,]$ref_name
    truth_mat <- cyto.vals$truth[[valType]][[valDataset]]
    truth_mat <- round(truth_mat, round_results)


    cyto.Res.tmp <- cyto.Res %>%
      filter(ref_name == refName & val_dataset == valDataset)

    yy <- cyto.Res.tmp[cyto.Res.tmp$method == "BayesPrism",]
    yy$method <- "xCell2"
    yy$res <- pull(filter(benchmark_table, ref_name == refName & val_dataset == valDataset), res)
    cyto.Res.tmp <- rbind(yy, cyto.Res.tmp)

    shared_celltypes <- Reduce(intersect, lapply(cyto.Res.tmp$res, rownames))
    shared_samples <- Reduce(intersect, lapply(cyto.Res.tmp$res, colnames))



    # ref = cyto.Res.tmp[1,]$ref_name
    # res = cyto.Res.tmp[1,]$res[[1]]
    # truth = truth_mat
    # shared_cts = shared_celltypes
    # shared_samp = shared_samples
    # get_spill_cors = spillcors
    # cor_method = cMethod


    out <- cyto.Res.tmp %>%
      rowwise() %>%
      mutate(cors = list(getCors(ref = ref_name, res, truth = truth_mat, shared_cts = shared_celltypes, shared_samp = shared_samples, get_spill_cors = spillcors, cor_method = cMethod))) %>%
      dplyr::select(method, cors) %>%
      unnest(cols = c(cors)) %>%
      mutate(ref = refName,
             val = valDataset)

    return(out)


  }, mc.cores = 15) %>%
    bind_rows()

  if (!weight_rhos) {
    return(all_cors)
  }

  if (by_val) {

    if (!is.null(ref2use)) {
      all_cors <- all_cors %>%
        filter(ref %in% ref2use)
    }

    all_cors_ref_combined <- all_cors %>%
      group_by(method, ref, val) %>%
      dplyr::summarise(cors_list = list(cor),
                       n_ct_samples = list(n)) %>% # Rhos are weighted by number of samples per cell type
      rowwise() %>%
      mutate(n_ct_samples = list(ifelse(n_ct_samples > 30, 30, n_ct_samples))) %>%
      mutate(ref_rho = list(combineRhos(rhos = cors_list, sample_sizes = log(n_ct_samples), use_median = FALSE, summarize = FALSE)),
             n_val_cts = length(cors_list))

    data_combined <- all_cors_ref_combined %>%
      mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
      unnest(ref_rho) %>%
      ungroup() %>%
      mutate(method = factor(method),
             val = factor(val)) %>%
      dplyr::select(-c(cors_list, n_ct_samples))

  }else{

    if (!is.null(ref2use)) {
      all_cors <- all_cors %>%
        filter(ref %in% ref2use)
    }

    all_cors_ref_combined <- all_cors %>%
      group_by(method, ref) %>%
      dplyr::summarise(cors_list = list(cor),
                       n_ct_samples = list(n)) %>% # Rhos are weighted by number of samples per cell type
      rowwise() %>%
      mutate(n_ct_samples = list(ifelse(n_ct_samples > 30, 30, n_ct_samples))) %>%
      mutate(ref_rho = list(combineRhos(rhos = cors_list, sample_sizes = log(n_ct_samples), use_median = FALSE, summarize = FALSE)),
             n_val_cts = length(cors_list))

    data_combined <- all_cors_ref_combined %>%
      mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
      unnest(ref_rho) %>%
      ungroup() %>%
      mutate(method = factor(method),
             ref = factor(ref)) %>%
      dplyr::select(-c(cors_list, n_ct_samples))

    return(data_combined)

  }





}





