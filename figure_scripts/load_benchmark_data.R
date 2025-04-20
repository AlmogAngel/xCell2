library(tidyverse)
library(xCell2)
library(BiocParallel)

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


methods2use <- c("xCell2", "BayesPrism", "CIBERSORTx", "EPIC", "MCPcounter", "dtangle", "Bisque", "DWLS", "SCDC", "DeconRNASeq", "Scaden")
mouse <- FALSE


# Mouse Benchmark -----

if (mouse) {
  # Load mouse benchmark data (prep_benchmark_data.R)
  mouse.refs <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/mouse_references.rds")
  mouse.vals <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/mouse_validation.rds")
  refs.vals.matched.mouse <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/refs.vals.matched.mouse.rds")
  
  # Load mouse cell type estimation results from other methods (run_benchmark_others.R)
  files <- list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/mouse", pattern = ".cyto.res.rds", full.names = TRUE)
  mouse.predicted <- lapply(files, function(f){
    
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
  mouse.predicted[mouse.predicted$method == "dwlr", ]$method <- "DWLR"
  mouse.predicted[mouse.predicted$method == "scdc", ]$method <- "SCDC"
  mouse.predicted[mouse.predicted$method == "bisque", ]$method <- "Bisque"
  
  # Set xCell2 for benchmark
  refs.vals.matched.mouse.xcell2 <- refs.vals.matched.mouse %>%
    mutate(n_val_samples = ncol(mouse.vals$truth[[val_type]][[val_dataset]])) %>%
    mutate(method = "xCell2", .before = everything())
}


# Human ------------

if (!mouse) {
  # Load human benchmark data (prep_benchmark_data.R)
  human.refs <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/human_references.rds")
  human.vals <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/human_validation.rds")
  refs.vals.matched.human <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/refs.vals.matched.human.rds")
  
  # Load human cell type estimation results from other methods (run_benchmark_others.R)
  files <- list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions", pattern = ".cyto.res.rds", full.names = TRUE)
  human.predicted <- lapply(files, function(f){
    
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
  human.predicted[human.predicted$method == "dwlr", ]$method <- "DWLR"
  human.predicted[human.predicted$method == "scdc", ]$method <- "SCDC"
  human.predicted[human.predicted$method == "bisque", ]$method <- "Bisque"
  
  # Set xCell2 for benchmark
  refs.vals.matched.human.xcell2 <- refs.vals.matched.human %>%
    mutate(n_val_samples = ncol(human.vals$truth[[val_type]][[val_dataset]])) %>%
    mutate(method = "xCell2", .before = everything())
}

# Benchmark functions ------------

vals.refs.res <- if(mouse) {refs.vals.matched.mouse.xcell2} else {refs.vals.matched.human.xcell2}
cyto.vals <- if(mouse) {mouse.vals} else {human.vals}
refsRDSList <- if(mouse) {mouse.refs} else {human.refs}
cyto.Res <- if(mouse) {mouse.predicted} else {human.predicted}

get_xcell2_benchmark_results <- function(cyto.vals = cyto.vals, benchmark_table = vals.refs.res, save_object = TRUE, dir = NULL, params, ncores, output_name){
  
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
      parallel_param <- MulticoreParam(workers = params$num_threads)
      xcell2_object <- xCell2Train(ref = ref.in, mix = mix.in, labels = labels, refType = refType,
                                   lineageFile = lineage_file, BPPARAM = parallel_param, 
                                   useOntology = params$use_ontolog, returnSignatures = params$return_signatures,
                                   returnAnalysis = params$return_analysis, useSpillover = params$use_sillover,
                                   spilloverAlpha = params$spillover_alpha, minPbCells = params$min_pb_cells,
                                   minPbSamples = params$min_pb_samples, minScGenes = params$min_sc_genes)
      saveRDS(xcell2_object, file)
    }
    
    if (!params$return_analysis) {
      parallel_param <- MulticoreParam(workers = params$num_threads)
      res <- xCell2Analysis(mix = mix.in, xcell2object = xcell2_object, rawScores = !params$use_sillover,
                            spillover = params$use_sillover, spilloverAlpha = params$spillover_alpha, BPPARAM = parallel_param)
    }
    
    
    return(res)
    
  }, mc.cores = ncores)
  saveRDS(xCell2results, output_name)
  
  print(paste0("Done - xCell2 benchmarking analysis results located in: ", output_name))
  
}


get_xcell2_correlations <- function(benchmark_table = vals.refs.res,  xCell2results, round_results = 3, weight_cors = TRUE,
                                    by_val = FALSE, ref2use = NULL, cMethod = "spearman", just_xcell2 = FALSE){
  
  getCors <- function(ref, res, truth, shared_cts, shared_samp, cor_method){
    
    # Get shared samples and cell types between predictions and ground truth
    celltypes <- intersect(rownames(res), rownames(truth))
    celltypes <- celltypes[celltypes %in% shared_cts]
    samples <- intersect(colnames(res), colnames(truth))
    samples <- samples[samples %in% shared_samp]
    
    # Generate a tibble with matched ground truth and prediction results for each sample
    df <- lapply(celltypes, function(ct){
      truth <- truth[ct,samples]
      res <- res[ct,samples]
      tibble(celltype = ct, truth = truth, prediction = res)
    }) %>%
      bind_rows()
    
    # Calculate correlations
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
  weightCors <- function(cors, sample_sizes = NULL, calc_median = FALSE){
    
    if (length(cors) == 1) {
      return(cors)
    }
    
    if (is.null(sample_sizes)) {
      sample_sizes <- rep(1, length(cors))
    }
    
    if (length(sample_sizes) != length(cors)) {
      stop("values and weights must have the same length")
    }
    
    # Fisher's Z Transformation
    cors[cors == 1] <- 0.999999999
    cors[cors == -1] <- -0.999999999
    z_values <- 0.5 * log((1 + cors) / (1 - cors))
    
    # Weight z-values
    z_weighted <- sample_sizes * z_values / mean(sample_sizes)
    
    if (calc_median) {
      z_weighted <- median(z_weighted)
    }
    
    # Inverse Fisher transformation
    cors_weighted <- (exp(2 * z_weighted) - 1) / (exp(2 * z_weighted) + 1)
    
    return(cors_weighted)
  }
  
  
  # Round data
  if (!just_xcell2) {
    cyto.Res <- cyto.Res %>%
      rowwise() %>%
      mutate(res = list(round(res, round_results)))
  }

  
  if (!is.null(round_results)) {
    xCell2results <- lapply(xCell2results, function(x){round(x, round_results)})
  }
  
  # Assign xCell2 results
  benchmark_table$res <- xCell2results
  
  if (just_xcell2) {
    ref_val_pairs <- benchmark_table |> 
      select(ref_tissue, ref_name, val_type, val_dataset)
    
    cyto.Res <- benchmark_table
  }else{
    # Calculate correlations for each pair of reference-validation data
    ref_val_pairs <- cyto.Res %>%
      group_by(ref_tissue, ref_name, val_type, val_dataset) %>%
      summarise(n = n()) %>%
      filter(n == length(methods2use)-1) %>%
      dplyr::select(-n) %>%
      ungroup()
    
  }
  

  
  
  all_cors <- parallel::mclapply(1:nrow(ref_val_pairs), function(i){
    
    # Get reference-validation data
    valType <- ref_val_pairs[i,]$val_type
    valDataset <- ref_val_pairs[i,]$val_dataset
    refName <- ref_val_pairs[i,]$ref_name
    
    # Load ground truth
    truth_mat <- cyto.vals$truth[[valType]][[valDataset]]
    truth_mat <- round(truth_mat, round_results)
    
    # Subset referemce-validation pair
    cyto.Res.tmp <- cyto.Res %>%
      filter(ref_name == refName & val_dataset == valDataset)
    
    # Add xCell2 results
    yy <- cyto.Res.tmp[1,]
    yy$method <- "xCell2"
    yy$res <- pull(filter(benchmark_table, ref_name == refName & val_dataset == valDataset), res)
    cyto.Res.tmp <- rbind(yy, cyto.Res.tmp)
    
    # Get shared samples and cell types
    shared_celltypes <- Reduce(intersect, lapply(cyto.Res.tmp$res, rownames))
    shared_samples <- Reduce(intersect, lapply(cyto.Res.tmp$res, colnames))
    
    # Calculate correlations for all methods
    out <- cyto.Res.tmp %>%
      rowwise() %>%
      mutate(cors = list(getCors(ref = ref_name, res, truth = truth_mat, shared_cts = shared_celltypes, shared_samp = shared_samples, cor_method = cMethod))) %>%
      dplyr::select(method, cors) %>%
      unnest(cols = c(cors)) %>%
      mutate(ref = refName,
             val = valDataset)
    
    return(out)
    
    
  }, mc.cores = 15) %>%
    bind_rows()
  
  if (!weight_cors) {
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
                       n_ct_samples = list(n)) %>% # Cors are weighted by number of samples per cell type
      rowwise() %>%
      mutate(n_ct_samples = list(ifelse(n_ct_samples > 30, 30, n_ct_samples))) %>% # Maximum weight is 30 samples
      mutate(ref_val_cor_median = weightCors(cors = cors_list, sample_sizes = log(n_ct_samples), calc_median = TRUE), # Natural logarithm of the number of samples
             n_val_cts = length(cors_list))
    
    data_combined <- all_cors_ref_combined %>%
      mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
      mutate(method = factor(method),
             val = factor(val)) %>%
      dplyr::select(-c(cors_list, n_ct_samples))
    
  }else{
    
    # Choose references to for weighted correlation calculation
    if (!is.null(ref2use)) {
      all_cors <- all_cors %>%
        filter(ref %in% ref2use)
    }
    
    # Calculate weighted correlation 
    all_cors_ref_combined <- all_cors %>%
      group_by(method, ref) %>%
      dplyr::summarise(cors_list = list(cor),
                       n_ct_samples = list(n)) %>% # Rhos are weighted by number of samples per cell type
      rowwise() %>%
      mutate(n_ct_samples = list(ifelse(n_ct_samples > 30, 30, n_ct_samples))) %>% # Maximum weight is 30 samples
      mutate(ref_cor = list(weightCors(cors = cors_list, sample_sizes = log(n_ct_samples))),  # Natural logarithm of the number of samples
             n_val_cts = length(cors_list))
    
    
    data_combined <- all_cors_ref_combined %>%
      mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
      unnest(ref_cor) %>%
      ungroup() %>%
      mutate(method = factor(method),
             ref = factor(ref)) %>%
      dplyr::select(-c(cors_list, n_ct_samples))
    
    return(data_combined)
    
  }
  
  
  
  
  
}





