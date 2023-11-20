library(tidyverse)
library(xCell2)
library(parallel)

# Load reference
ref.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/lm22_ref.rds")
ref = ref.in$ref
labels = ref.in$labels

# Load mixture
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
val_dataset = "BG_blood"
val_type = "blood"
mix <- cyto.vals$mixtures[[val_type]][[val_dataset]]

# Set data type
data_type = "array"

# Get shared genes
if (data_type == "sc") {
  shared_clean_genes <- xCell2CleanGenes(ref = ref, mix = mix, top_var_genes = TRUE, use_protein_coding = TRUE, n_var_genes = 5000)
}else{
  shared_clean_genes <- xCell2CleanGenes(ref = ref, mix = mix, top_var_genes = FALSE, use_protein_coding = FALSE, n_var_genes = 5000)
}
ref <- shared_clean_genes$ref
mix <- shared_clean_genes$mix

# Load parameters
lineage_file = ref.in$lineage_file
sim_fracs = c(0, seq(0.01, 0.25, 0.002), seq(0.3, 1, 0.05))
diff_vals = c(1, 1.585, 2, 3, 4, 5)
probs = c(0.01, 0.05, 0.1, 0.25, 0.333, 0.49)
min_genes = 3
max_genes = 100
return_sigs = FALSE
sigsFile = NULL
minPBcells = 30
minPBsamples = 10
weightGenes = TRUE
medianGEP = TRUE
sim_noise = NULL
ct_sims = 20
sims_sample_frac = 0.1
seed = 123
nCores = 10
# mix = NULL
simMethod = "ref_mix_thin"
# c("ref_multi", "ref_thin", "ref_mix_thin")


# Load signature
# sigsFile = "/bigdata/almogangel/xCell2_data/dev_data/sigs/BG_blood_ts_blood_sigs.rds"


validateInputs <- function(ref, labels, data_type){
  if (length(unique(labels$label)) < 3) {
    stop("Reference must have at least 3 cell types")
  }

  if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
    stop("ref must be one of those classes: matrix, dgCMatrix, Matrix")
  }

  if (!"data.frame" %in% class(labels)) {
    stop("labels must be a dataframe.")
  }

  if (!data_type %in% c("rnaseq", "array", "sc")) {
    stop("data_type should be 'rnaseq', 'array' or 'sc'.")
  }

  if (sum(grepl("_", labels$label)) != 0) {
    message("Changing underscores to dashes in cell-types labels!")
    labels$label <- gsub("_", "-", labels$label)
  }


  out <- list(ref = ref,
              labels = labels)
  return(out)

}
# TODO: Write a function for TPM/CPM normalization for bulk reference
normRef <- function(ref, data_type){

  if(all(ref == floor(ref))){
    if (data_type == "sc") {
      message("Normalizing scRNA-Seq counts to CPM.")
      lib_sizes <- Matrix::colSums(ref)
      norm_factor <- 1000000 / lib_sizes
      ref_norm <- ref %*% Matrix::Diagonal(x = norm_factor)
      colnames(ref_norm) <- colnames(ref)
      ref_norm <- as.matrix(ref_norm)
      return(ref_norm)
    }else{
      # TODO: Write a function for TPM normalization for bulk reference
      message("Normalizing counts to TPM.")
      return(ref_norm)
    }
  }else{
    message("Assuming reference already normalized.")
    return(ref)
  }

}
sc2pseudoBulk <- function(ref, labels, min_n_cells, min_ps_samples, seed2use){

  set.seed(seed2use)

  celltypes <- unique(labels$label)

  groups_list <- lapply(celltypes, function(ctoi){

    ctoi_samples <- labels[labels$label == ctoi,]$sample

    # Calculate maximum possible number of groups given min_n_cells
    num_groups <- ceiling(length(ctoi_samples) / min_n_cells)
    if (num_groups < min_ps_samples) {
      num_groups <- min_ps_samples
    }

    # Generate min_ps_samples pseudo samples of CTOI
    if (length(ctoi_samples) > min_ps_samples) {

      ctoi_samples_shuffled <- sample(ctoi_samples, length(ctoi_samples))
      list_of_ctoi_samples_shuffled <- split(ctoi_samples_shuffled, ceiling(seq_along(ctoi_samples_shuffled) / (length(ctoi_samples_shuffled) / num_groups)))

      sapply(list_of_ctoi_samples_shuffled, function(ctoi_group){
        if (length(ctoi_group) == 1) {
          ref[,ctoi_group]
        }else{
          if("matrix" %in% class(ref)) Rfast::rowsums(ref[,ctoi_group]) else Matrix::rowSums(ref[,ctoi_group])
        }
      })

    }else{
      ref[,ctoi_samples]
    }

  })

  names(groups_list) <- celltypes
  pseudo_ref <- as.matrix(bind_cols(groups_list))
  rownames(pseudo_ref) <- rownames(ref)

  pseudo_label <- tibble(labels) %>%
    select(ont, label) %>%
    unique() %>%
    right_join(., tibble(label = sub("\\.\\d+$", "", colnames(pseudo_ref)), sample = colnames(pseudo_ref), dataset = "pseudoBulk"), by = "label") %>%
    as.data.frame()

  return(list(ref = pseudo_ref, labels = pseudo_label))

}
makeGEPMat <- function(ref, labels, use_median){

  celltypes <- unique(labels$label)

  gep_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      if(use_median){
        type_vec <- Rfast::rowMedians(as.matrix(ref[,type_samples]))
      }else{
        type_vec <- Rfast::rowmeans(as.matrix(ref[,type_samples]))
      }
    }
  })
  rownames(gep_mat) <- rownames(ref)

  return(gep_mat)
}
getCellTypeCorrelation <- function(gep_mat, data_type){

  celltypes <- colnames(gep_mat)

  if (data_type != "sc") {

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
logTransformRef <- function(ref){

  if(max(ref) >= 50){
    message("> Transforming reference to log2-space (maximum expression value >= 50).")
    ref.log2 <- log2(ref+1)
    return(ref.log2)
  }else{
    message("> Assuming reference is already in log2-space (maximum expression value < 50).")
    return(ref)
  }

}
makeQuantiles <- function(ref, labels, probs, ncores){

  celltypes <- unique(labels[,2])

  quantiles_mat_list <-  parallel::mclapply(celltypes, function(type){

    type_samples <- labels[,2] == type

    if (sum(type_samples) == 1) {
      # If there is one sample for this cell type -> duplicate the sample to make a data frame
      type.df <- cbind(ref[,type_samples], ref[,type_samples])
    }else{
      type.df <- ref[,type_samples]
    }

    # Calculate quantiles
    # TODO: Balance quantiles by dataset
    type_quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1-probs))), na.rm=TRUE))
  }, mc.cores = ncores, mc.set.seed = FALSE)
  names(quantiles_mat_list) <- celltypes

  return(quantiles_mat_list)
}
createSignatures <- function(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes, ncores){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Set signature thresholds grid
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)
    param.df <- param.df[order(-param.df$diff_vals, param.df$probs), ]

    # Generate signatures
    type_sigs <- list()
    for (i in 1:nrow(param.df)){


      # Get a Boolean matrices with genes that pass the quantiles criteria
      diff <- param.df[i, ]$diff_vals # difference threshold
      lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile cutoff
      upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile cutoff

      diff_genes.mat <- sapply(not_dep_celltypes, function(x){
        get(type, quantiles_matrix)[lower_prob,] > get(x, quantiles_matrix)[upper_prob,] + diff
      })


      if(weight_genes){
        # Score genes using cell types correlations as weights
        type_weights <- cor_mat[type, not_dep_celltypes]
        type_weights[type_weights < 0.001] <- 0.001 # Fix minimum correlation to avoid zero and negative correlations
        gene_scores <- apply(diff_genes.mat, 1, function(x){
          sum(type_weights[which(x)])
        })
      }else{
        # All cell types scores are the same
        gene_scores <- apply(diff_genes.mat, 1, function(x){
          sum(x)
        })
      }

      gene_passed <- sort(gene_scores[gene_scores > 0], decreasing = TRUE)

      # If less than min_genes passed move to next parameters
      if (length(gene_passed) < min_genes) {
        next
      }

      # Round and sort top genes scores
      top_scores <- sort(unique(round(gene_passed-0.5)), decreasing = TRUE)

      # Take top 3 highest scores from top_scores
      top_scores <- top_scores[1:3]

      for (score in top_scores) {

        n_genes <- sum(gene_passed >= score)

        if (n_genes < min_genes) {
          next
        }

        if (n_genes > max_genes) {
          break
        }

        sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, sep = "_")
        type_sigs[[sig_name]] <- names(which(gene_passed >= score))
      }

    }


    # Remove duplicate signatures
    type_sigs_sorted <- lapply(type_sigs, function(x) sort(x))
    type_sigs_sorted_collapsed <- sapply(type_sigs_sorted, paste, collapse = ",")
    duplicated_sigs <- duplicated(type_sigs_sorted_collapsed)
    type_sigs <- type_sigs[!duplicated_sigs]

    if (length(type_sigs) < 3) {
      warnings(paste0("Not enough signatures found for ", type))
    }

    return(type_sigs)
  }


  celltypes <- unique(labels[,2])

  all_sigs <- parallel::mclapply(celltypes, function(type){
    getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes)
  }, mc.cores = ncores, mc.set.seed = FALSE)


  all_sigs <- unlist(all_sigs, recursive = FALSE)


  if (length(all_sigs) == 0) {
    stop("No signatures found for reference!")
  }


  return(all_sigs)
}
addEssentialGenes <- function(ref, signatures){

  data("celltype.data", package = "xCell2")
  celltypes <- unique(gsub("#.*", "", names(signatures)))

  for (ct in celltypes) {

    # Check if cell type exists in celltype.data
    if (!ct %in% celltype.data$all_labels) {
      next
    }

    # Get essential genes
    ct_label <- celltype.data[celltype.data$all_labels == ct,]$xCell2_labels
    essen_genes <- unique(unlist(celltype.data[celltype.data$xCell2_labels == ct_label,]$essential_genes))

    if (all(is.na(essen_genes))) {
      next
    }

    essen_genes <- intersect(rownames(ref), essen_genes) # Use only essential genes which are in ref

    # Add to signature
    ct_sigs <- which(startsWith(names(signatures), paste0(ct, "#")))
    for (sig in ct_sigs) {
      signatures[sig][[1]] <- unique(c(signatures[sig][[1]], essen_genes))
    }

  }

  return(signatures)

}
makeSimulations <- function(ref, labels, mix, gep_mat, cor_mat, dep_list, sim_fracs, sim_method, ctoi_samples_frac, n_sims, ncores, seed2use){

  set.seed(seed2use)

  celltypes <- unique(labels$label)

  getSubMatrix <- function(mat, sim_fracs, n_samples_sim){
    if (class(mat)[1] == "numeric") {
      mat_sub <- matrix(rep(mat, length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs))
      rownames(mat_sub) <- rownames(mat)
    }else{
      mat_sub <- sapply(1:length(sim_fracs), function(i){
        mat_tmp <- mat[,sample(1:ncol(mat), n_samples_sim)]
        if (class(mat_tmp)[1] == "numeric") {
          mat_tmp
        }else{
          rowMeans(mat_tmp)
        }
      })
    }

    rownames(mat_sub) <- rownames(mat)

    return(mat_sub)
  }
  adjustLibSize <- function(ctoi_mat, controls_mat){

    # Scale to simulate counts data
    scale_factor <- 10000
    ctoi_mat_scaled <- round(ctoi_mat * scale_factor)
    controls_mat_scaled <- round(controls_mat * scale_factor)

    # Adjust reference-controls library size
    min_lib_size <- min(min(colSums(ctoi_mat_scaled)), min(colSums(controls_mat_scaled)))
    ref_ctoi_sub_lib_fracs <- min_lib_size/colSums(ctoi_mat_scaled)
    ref_controls_sub_lib_fracs <- min_lib_size/colSums(controls_mat_scaled)

    # Thin data to adjust lib size
    ctoi_mat_scaled_thin <- seqgendiff::thin_lib(ctoi_mat_scaled, thinlog2 = -log2(ref_ctoi_sub_lib_fracs), type = "thin")$mat
    controls_mat_scaled_thin <- seqgendiff::thin_lib(controls_mat_scaled, thinlog2 = -log2(ref_controls_sub_lib_fracs), type = "thin")$mat

    # Unscale
    ctoi_mat_thin <- ctoi_mat_scaled_thin/scale_factor
    rownames(ctoi_mat_thin) <- rownames(ctoi_mat)
    controls_mat_thin <- controls_mat_scaled_thin/scale_factor
    rownames(controls_mat_thin) <- rownames(controls_mat)

    return(list(ctoi_mat_thin = ctoi_mat_thin, controls_mat_thin = controls_mat_thin))
  }
  makeFractionMatrix <- function(mat, sim_fracs, sim_method, control){


    # Check if data not in counts (all integers) because you can't thin fractions (?)
    if (sim_method != "ref_multi") {
      scale_factor <- 10000
      mat <- round(mat * scale_factor)
    }

    # Adjust simulation fractions for controls
    if (control) {
      sim_fracs <- 1-sim_fracs
    }


    # Multiply reference matrix to simulate fractions
    if (sim_method == "ref_multi") {
      sim <- mat %*% diag(sim_fracs)
    }else{
      # Thin reference/mixture to simulate fractions
      if (sum(sim_fracs == 0) != 0) { # Can't thin when frac = 0
        zero_index <- which(sim_fracs == 0)
        sim_fracs[zero_index] <- 0.001
        sim <- seqgendiff::thin_lib(mat, thinlog2 = -log2(sim_fracs), type = "thin")$mat
        sim_fracs[zero_index] <- 0
        sim <- sim/scale_factor
        sim[,zero_index] <- 0
      }else{
        sim <- seqgendiff::thin_lib(mat, thinlog2 = -log2(sim_fracs), type = "thin")$mat
      }
    }

    rownames(sim) <- rownames(mat)

    return(sim)
  }


  sim_list <- parallel::mclapply(celltypes, function(ctoi){

    ref_ctoi <- ref[,labels$label == ctoi]

    # Number of CTOI samples to use for each simulation
    n_samples_sim <- round(ncol(ref_ctoi) * ctoi_samples_frac)
    n_samples_sim <- ifelse(n_samples_sim < 1, 1, n_samples_sim)

    if (sim_method != "ref_mix_thin") {
      # Get control cell types
      dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
      controls <- celltypes[!celltypes %in% dep_cts]
    }

    # Generate n_sims simulations
    ctoi_sim_list <- lapply(1:n_sims, function(i){

      # Use n_samples_sim random samples for CTOI
      ref_ctoi_sub <- getSubMatrix(mat = ref_ctoi, sim_fracs, n_samples_sim)

      # Use n_samples_sim random samples for controls
      if (sim_method != "ref_mix_thin") {

        controls2use <- sample(controls, sample(1:length(controls), 1), replace = FALSE)
        samples2use <- labels %>%
          filter(label %in% controls2use) %>%
          pull(sample)
        ref_controls <- ref[,samples2use]

        n_samples_sim <- round(ncol(ref_controls) * ctoi_samples_frac)
        n_samples_sim <- ifelse(n_samples_sim < 1, 1, n_samples_sim)

        ref_controls_sub <- getSubMatrix(mat = ref_controls, sim_fracs, n_samples_sim)
      }else{

        # Shuffle expression values between genes
        mix_shuffled <- t(apply(mix, 1, sample))
        mix_shuffled <- mix_shuffled[rownames(mix),]

        # Number of control samples to use
        n_samples_sim <- round(ncol(mix_shuffled) * ctoi_samples_frac)
        n_samples_sim <- ifelse(n_samples_sim < 1, 1, n_samples_sim)

        ref_controls_sub <- getSubMatrix(mat = mix_shuffled, sim_fracs, n_samples_sim)
      }

      # Thin to adjust library size
      data_adjusted <- adjustLibSize(ctoi_mat = ref_ctoi_sub, controls_mat = ref_controls_sub)
      ref_ctoi_sub <- data_adjusted$ctoi_mat_thin
      ref_controls_sub <- data_adjusted$controls_mat_thin

      # Make fraction matrices
      ctoi_frac_mat <- makeFractionMatrix(mat = ref_ctoi_sub, sim_fracs, sim_method, control = FALSE)
      control_frac_mat <- makeFractionMatrix(mat = ref_controls_sub, sim_fracs, sim_method, control = TRUE)

      # Combine CTOI and control(s) fractions matrix
      simulation <- ctoi_frac_mat + control_frac_mat
      colnames(simulation) <- paste0("mix", "%%", sim_fracs)

      simulation


    })



    do.call(cbind, ctoi_sim_list)


  }, mc.cores = ncores, mc.set.seed = FALSE)
  names(sim_list) <- celltypes

  return(sim_list)

}
scoreSimulations <- function(signatures, simulations, ncores){

  celltypes <- names(simulations)

  sims_scored <- parallel::mclapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]
    ctoi_sim <- simulations[[ctoi]]
    sim_ranked <- singscore::rankGenes(ctoi_sim)
    colnames(sim_ranked) <- make.unique(colnames(sim_ranked), sep = "%")

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })

    score <- cbind(scores, frac = as.numeric(gsub("mix%%", "", colnames(ctoi_sim))))
    score

  }, mc.cores = ncores, mc.set.seed = FALSE)

  names(sims_scored) <- celltypes

  return(sims_scored)

}
trainModels <- function(simulations_scored, ncores, seed2use){

  set.seed(seed2use)

  fitModel <- function(data, nRFcores){

    options(rf.cores=nRFcores, mc.cores=1)
    model <- randomForestSRC::var.select(frac ~ ., as.data.frame(data), method = "vh.vimp", verbose = FALSE, refit = TRUE, fast = TRUE)
    selected_features <- model$topvars

    return(tibble(model = list(model$rfsrc.refit.obj), sigs_filtered = list(selected_features)))

  }

  if (ncores < 4) {
    mcCores <- round(ncores*(3/4))
    rfCores <- round(ncores*(1/4))
  }else{
    mcCores <- ncores
    rfCores <- 1
  }

  #start <- Sys.time()
  models_list <- parallel::mclapply(simulations_scored, function(data){
    fitModel(data, rfCores)
  }, mc.cores = mcCores, mc.set.seed = FALSE)
  #end <- Sys.time()
  #print(end-start)


  enframe(models_list, name = "celltype") %>%
    unnest(value) %>%
    return(.)

}
getSpillOverMat <- function(simulations, signatures, dep_list, models, frac2use){

  scoreTransform <- function(mat, signatures, models){

    # Score
    mat_ranked <- singscore::rankGenes(mat)
    scores <- sapply(signatures, simplify = TRUE, function(sig){
      singscore::simpleScore(mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    rownames(scores) <- colnames(mat)


    transfomed_tbl <- models %>%
      rowwise() %>%
      mutate(predictions = list(round(randomForestSRC::predict.rfsrc(model, newdata = as.data.frame(scores[,startsWith(colnames(scores), paste0(celltype, "#"))]))$predicted, 4))) %>%
      select(celltype, predictions) %>%
      unnest_longer(predictions, indices_to = "sim_celltype") %>%
      pivot_wider(names_from = sim_celltype, values_from = predictions)


    # Convert to matrix
    transfomed_mat <- as.matrix(transfomed_tbl[,-1])
    rownames(transfomed_mat) <- pull(transfomed_tbl[,1])
    colnames(transfomed_mat) <- rownames(transfomed_mat)


    return(transfomed_mat)

  }


  # Get CTOIs matrix with frac2use fraction
  frac_col <- which(endsWith(colnames(simulations[[1]]), paste0("%%", frac2use)))
  ctoi_mat <- sapply(simulations, function(sim){
    if (length(frac_col) > 1) {
      apply(sim[,frac_col], 1, median)
    }else{
      sim[,frac_col]
    }
  })

  # Get control matrix with CTOI fraction = 0
  frac_col <- which(endsWith(colnames(simulations[[1]]), paste0("%%", 0)))
  controls_mat <- sapply(simulations, function(sim){
    if (length(frac_col) > 1) {
      apply(sim[,frac_col], 1, median)
    }else{
      sim[,frac_col]
    }
  })

  # Score and transform simulations
  sim_transformed <- scoreTransform(mat = ctoi_mat, signatures, models)
  controls_mat_transformed <- scoreTransform(mat = controls_mat, signatures, models)

  # Remove control signal from the transformed mixture
  spill_mat <- sim_transformed - controls_mat_transformed

  # Clean and normalize spill matrix
  spill_mat[spill_mat < 0] <- 0
  spill_mat <- spill_mat / diag(spill_mat)

  # Insert zero to dependent cell types
  for(ctoi in rownames(spill_mat)){
    dep_cts <- unname(unlist(dep_list[[ctoi]]))
    dep_cts <- dep_cts[dep_cts != ctoi]
    spill_mat[ctoi, dep_cts] <- 0
  }

  # TODO: Check this parameter
  spill_mat[spill_mat > 0.5] <- 0.5
  diag(spill_mat) <- 1

  # pheatmap::pheatmap(spill_mat, cluster_rows = F, cluster_cols = F)

  return(spill_mat)

}


#' @slot signatures list of xCell2 signatures
#' @slot dependencies list of cell type dependencies
#' @slot models data frame of cell type transformation models
#' @slot spill_mat matrix of cell types spillover
#' @slot genes_used character vector of genes names used to train the signatures
#' @importFrom methods new
# Create S4 object for the new reference
setClass("xCell2Signatures", slots = list(
  signatures = "list",
  dependencies = "list",
  models = "data.frame",
  spill_mat = "matrix",
  genes_used = "character"
))


#' xCell2Train function
#'
#' This function generates signatures for each cell type.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import readr
#' @import parallel
#' @importFrom randomForestSRC var.select
#' @importFrom Rfast rowMedians rowmeans rowsums
#' @importFrom parallel mclapply
#' @importFrom pbapply pblapply pbsapply
#' @importFrom RRF RRF
#' @importFrom seqgendiff thin_lib
#' @importFrom Matrix rowMeans rowSums colSums
#' @importFrom singscore rankGenes simpleScore
#' @param ref A reference gene expression matrix.
#' @param labels A data frame in which the rows correspond to samples in the ref. The data frame must have four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545" or NA if there is no ontology).
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample should match the column name in ref.
#'   "dataset": the cell type sample dataset or subject (for single-cell) as a character.
#' @param data_type Gene expression data type: "rnaseq", "array", or "sc".
#' @param lineage_file Path to the cell type lineage file generated with `xCell2GetLineage` function (optional).
#' @param weightGenes description
#' @param sim_fracs A vector of mixture fractions to be used in signature filtering (optional).
#' @param probs A vector of probability thresholds to be used for generating signatures (optional).
#' @param diff_vals A vector of delta values to be used for generating signatures (optional).
#' @param min_genes The minimum number of genes to include in the signature (optional).
#' @param max_genes The maximum number of genes to include in the signature (optional).
#' @param sigsFile description
#' @param return_sigs description
#' @param minPBcells description
#' @param minPBgroups description
#' @param ct_sims description
#' @param sims_sample_frac description
#' @param nCores description
#' @param mix description
#' @param simMethod description
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, data_type, mix = NULL, lineage_file = NULL, weightGenes = TRUE, medianGEP = TRUE, seed = 123, probs = c(0.01, 0.05, 0.1, 0.25, 0.333, 0.49),
                        sim_fracs = c(0, seq(0.01, 0.25, 0.005), seq(0.3, 1, 0.05)), diff_vals = c(1, 1.585, 2, 3, 4, 5),
                        min_genes = 3, max_genes = 100, return_sigs = FALSE, sigsFile = NULL, minPBcells = 30, minPBsamples = 10,
                        ct_sims = 20, sims_sample_frac = 0.1, simMethod = "ref_thin", nCores = 1){


  # Validate inputs
  inputs_validated <- validateInputs(ref, labels, data_type)
  ref <- inputs_validated$ref
  labels <- inputs_validated$labels


  # TODO: first sum counts and then normalize or vice versa?

  # Generate pseudo bulk from scRNA-Seq reference
  if (data_type == "sc") {
    message("Converting scRNA-seq reference to pseudo bulk...")
    ps_data <- sc2pseudoBulk(ref, labels, min_n_cells = minPBcells, min_ps_samples = minPBsamples, seed2use = seed)
    ref <- ps_data$ref
    labels <- ps_data$labels
  }

  # Normalize reference
  ref <- normRef(ref, data_type)


  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  gep_mat <- makeGEPMat(ref, labels, use_median = medianGEP)
  cor_mat <- getCellTypeCorrelation(gep_mat, data_type)


  # Get cell type dependencies list
  message("Loading dependencies...")
  if (is.null(lineage_file)) {
    dep_list <- xCell2::xCell2GetLineage(labels, out_file = NULL)
  }else{
    dep_list <- getDependencies(lineage_file)
  }


  # Generate/Load signatures
  if (is.null(sigsFile)) {
    # Generate signatures
    message("Calculating quantiles...")
    ref_log <- logTransformRef(ref) # Log2-transformation
    quantiles_matrix <- makeQuantiles(ref_log, labels, probs, ncores = nCores)
    message("Generating signatures...")
    signatures <- createSignatures(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes = weightGenes, ncores = nCores)

    # Add essential genes
    signatures <- addEssentialGenes(ref, signatures)

    if (return_sigs) {
      return(signatures)
    }

  }else{
    # Load signatures
    message("Loading signatures...")
    signatures <- readRDS(sigsFile)
  }


  # Make simulations
  message("Generating simulations...")
  simulations <- makeSimulations(ref, labels, mix, gep_mat, cor_mat, dep_list, sim_fracs, sim_method = simMethod, ctoi_samples_frac = sims_sample_frac, n_sims = ct_sims, ncores = nCores, seed2use = seed)
  message("Scoring simulations...")
  simulations_scored <- scoreSimulations(signatures, simulations, nCores)


  # Filter signatures and train RF model
  message("Filtering signatures and training models...")
  models <- trainModels(simulations_scored, ncores = nCores, seed2use = seed)
  signatures <- signatures[unlist(models$sigs_filtered)]
  models <- models[,-3]


  # Get spillover matrix
  message("Generating spillover matrix...")
  frac2use <- sim_fracs[which.min(abs(sim_fracs - 0.25))]
  spill_mat <- getSpillOverMat(simulations, signatures, dep_list, models, frac2use)


  # Save results in S4 object
  xCell2Sigs.S4 <- new("xCell2Signatures",
                       signatures = signatures,
                       dependencies = dep_list,
                       models = models,
                       spill_mat = spill_mat,
                       genes_used = rownames(ref))


  return(xCell2Sigs.S4)

}
