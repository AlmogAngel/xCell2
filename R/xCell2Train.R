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
makeQuantiles <- function(ref, labels, probs, dep_list){

  celltypes <- unique(labels[,2])

  quantiles_mat_list <-  pbapply::pblapply(celltypes, function(type){

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
  })
  names(quantiles_mat_list) <- celltypes

  return(quantiles_mat_list)
}
createSignatures <- function(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes, ncores){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Set signature thresholds grid
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)

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

      # Take top n_top_scores highest scores
      # n_top_scores <- ifelse(length(top_scores) > n_top_scores, n_top_scores, length(top_scores))
      # for (score in top_scores[1:n_top_scores])

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
makeSimulations <- function(ref, labels, mix, gep_mat, cor_mat, dep_list, sim_fracs, sim_method, ctoi_samples_frac, n_sims, noise, ncores, seed2use){

  set.seed(seed2use)

  celltypes <- unique(labels$label)

  makeSamplesPool <- function(labels, ctoi){

    # Mix CTOI samples by datasets
    ctoi_samples_pool <- c()
    while(!all(labels[labels$label == ctoi,]$sample %in% ctoi_samples_pool)) {
      ctoi_samples_pool <- c(ctoi_samples_pool,
                             labels %>%
                               filter(label == ctoi & !sample %in% ctoi_samples_pool) %>%
                               slice_head(n = 1, by = dataset) %>%
                               pull(sample))
    }

    return(ctoi_samples_pool)
  }
  makeFractionMatrix <- function(expr_mat, sim_fracs, sim_method, control){

    # Get mean expression vector
    if (class(expr_mat)[1] == "numeric") {
      mean_expression <- expr_mat
    }else{
      mean_expression <- rowMeans(expr_mat)
    }

    # Check if data not in counts (all integers) because you can't thin fractions (?)
    if (sim_method != "ref_multi") {
      scale_factor <- 10000 # TODO: Ask Anna
      mean_expression <- round(mean_expression * scale_factor)
    }

    # Adjust simulation fractions for controls
    if (control) {
      sim_fracs <- 1-sim_fracs
    }

    # Build simulation matrix
    mat <- matrix(rep(mean_expression, length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs))
    rownames(mat) <- rownames(ref)


    # Multiply reference matrix to simulate fractions
    if (sim_method == "ref_multi") {
      sim <- mat %*% diag(sim_fracs)
    }else{
      # Thin reference/mixture to simulate fractions
      if (sum(sim_fracs == 0) != 0) { # Can't thin when frac = 0
        zero_index <- which(sim_fracs == 0)
        sim <- seqgendiff::thin_lib(mat[,-zero_index], thinlog2 = -log2(sim_fracs[-zero_index]), type = "thin")$mat
        sim <- sim/scale_factor
        if (zero_index == 1) {
          # Zero is  the first column
          sim <- cbind(mat[,1]*0, sim)
        }else{
          # Zero is the last column
          sim <- cbind(sim, mat[,1]*0)
        }
      }else{
        sim <- seqgendiff::thin_lib(mat, thinlog2 = -log2(sim_fracs), type = "thin")$mat
      }
    }

    return(sim)
  }


  sim_list <- parallel::mclapply(celltypes, function(ctoi){

    # Generate CTOI samples pool
    ctoi_samples_pool <- makeSamplesPool(labels, ctoi)

    # Number of CTOI samples to use
    n_samples_sim <- round(length(ctoi_samples_pool) * ctoi_samples_frac)
    n_samples_sim <- ifelse(n_samples_sim < 1, 1, n_samples_sim)

    if (sim_method != "ref_mix_thin") {
      # Get control cell types
      dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
      controls <- celltypes[!celltypes %in% dep_cts]
    }

    # Generate n_sims simulations
    ctoi_sim_list <- lapply(1:n_sims, function(i){
      # Make CTOI fraction matrix
      samples2use <- ctoi_samples_pool[1:n_samples_sim]
      ref_sub <- ref[,samples2use]
      ctoi_frac_mat <- makeFractionMatrix(expr_mat = ref_sub, sim_fracs, sim_method, control = FALSE)

      # Move samples2use to the end of the vector
      ctoi_samples_pool <- c(ctoi_samples_pool[!ctoi_samples_pool %in% samples2use], samples2use)

      # Make control(s) fractions matrix
      if (sim_method != "ref_mix_thin") {
        controls2use <- sample(controls, sample(1:length(controls), 1), replace = FALSE)
        samples2use <- labels %>%
          filter(label %in% controls2use) %>%
          group_by(label) %>%
          sample_n(1) %>%
          pull(sample)
        ref_sub <- ref[,samples2use]
        control_frac_mat <- makeFractionMatrix(expr_mat = ref_sub, sim_fracs, sim_method, control = TRUE)
      }else{
        # Shuffle expression values between genes
        mix_shuffled <- t(apply(mix, 1, sample))
        mix_shuffled <- mix_shuffled[rownames(mix),]
        # Select random samples
        num_columns_to_sample <- sample(2:ncol(mix_shuffled), 1)
        mix_shuffled <- mix_shuffled[, sample(1:ncol(mix_shuffled), num_columns_to_sample)]
        control_frac_mat <- makeFractionMatrix(mix_shuffled, sim_fracs, sim_method, control = TRUE)
      }


      # Combine CTOI and control(s) fractions matrix
      simulation <- ctoi_frac_mat + control_frac_mat
      colnames(simulation) <- paste0("mix", "%%", sim_fracs)

      simulation


    })

    # TODO: Noise
    if (!is.null(noise)) {
      ctoi_sim_list_noised <- lapply(ctoi_sim_list, function(sim){

        random_ref_samples <- sample(1:ncol(ref), size = ncol(sim), replace = FALSE)

        if (sum(sim_fracs == 0) != 0) {
          zero_index <- which(sim_fracs == 0)
          noiseLevel <- -log2(sim_fracs[-zero_index]/(1/noise)) # % noise is from the original frac
          ref_noise <- seqgendiff::thin_lib(round(ref[,random_ref_samples])[,-zero_index], thinlog2 = noiseLevel, type = "thin")$mat
          if (zero_index == 1) {
            # Zero is  the first column
            ref_noise <- cbind(sim[,1]*0, ref_noise)
          }else{
            # Zero is the last column
            ref_noise <- cbind(ref_noise, sim[,1]*0)
          }
        }else{
          ref_noise <- seqgendiff::thin_lib(round(ref[,random_ref_samples])[,-zero_index], thinlog2 = -log2(sim_fracs/(1/noise)), type = "thin")$mat
        }

        sim + ref_noise
      })
      do.call(cbind, ctoi_sim_list_noised)
    }else{
      do.call(cbind, ctoi_sim_list)
    }

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
trainModels <- function(simulations_scored, regGamma, ncores, seed2use){

  set.seed(seed2use)

  fitModel <- function(data, gamma){

    # Feature selection via GRRF
    # https://sites.google.com/site/houtaodeng/rrf?authuser=0
    RF <- RRF::RRF(x = data[,-ncol(data)], y = data[,ncol(data)], flagReg = 0, importance = TRUE, ntree = 1000)
    RF_imp <- RF$importance[,"%IncMSE"] / max(RF$importance[,"%IncMSE"])
    RRF <- RRF::RRF(x = data[,-ncol(data)], y = data[,ncol(data)], flagReg = 1, ntree = 1000, coefReg = (1-gamma) + gamma*RF_imp)
    selected_features <- colnames(data)[RRF$feaSet]
    data <- data[, c(selected_features, "frac")]


    # Build model
    # model <- RRF::tuneRRF(x = data[,-ncol(data)], y = data[,ncol(data)], flagReg = 0, importance = FALSE, ntreeTry=1000, stepFactor=1.1, improve=0.001, trace=FALSE, plot=FALSE, doBest=TRUE)
    model <- RRF::RRF(x = data[,-ncol(data)], y = data[,ncol(data)], flagReg = 0, ntree = 1000)

    return(tibble(model = list(model), sigs_filtered = list(selected_features)))
  }

  if (ncores == 1) {
    message("Using single core to train RF models (for faster run time use more cores).")
    models_list <- pbapply::pblapply(simulations_scored, function(data){
      fitModel(data, gamma = regGamma)
    })
  }else{
    models_list <- parallel::mclapply(simulations_scored, function(data){
      fitModel(data, gamma = regGamma)
    }, mc.cores = ncores, mc.set.seed = FALSE)
  }

  enframe(models_list, name = "celltype") %>%
    unnest(value) %>%
    return(.)

}
getSpillOverMat <- function(simulations, signatures, dep_list, trans_models, n_sims, frac2use){

  scoreTransform <- function(mat, signatures, trans_models, is_controls){

    # Score
    mat_ranked <- singscore::rankGenes(mat)
    scores <- sapply(signatures, simplify = TRUE, function(sig){
      singscore::simpleScore(mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    rownames(scores) <- colnames(mat)


    transfomed_tbl <- trans_models %>%
      mutate(predictions = list(round(predict(model, newdata = scores[,startsWith(colnames(scores), paste0(celltype, "#"))], type = "response"), 4))) %>%
      select(celltype, predictions) %>%
      unnest_longer(predictions, indices_to = "sim_celltype") %>%
      pivot_wider(names_from = sim_celltype, values_from = predictions)


    # Convert to matrix
    transfomed_mat <- as.matrix(transfomed_tbl[,-1])
    rownames(transfomed_mat) <- pull(transfomed_tbl[,1])

    if (!is_controls) {
      transfomed_mat <- transfomed_mat[colnames(mat), colnames(mat)]
    }

    return(transfomed_mat)

  }

  if (n_sims == 1) {

    # Get CTOIs  matrix with frac2use fraction
    frac_col <- which(endsWith(colnames(simulations[[1]]), paste0("%%", frac2use)))
    ctoi_mat <- sapply(simulations, function(sim){
      sim[,frac_col]
    })


    # Get control matrix with CTOI fraction = 0
    frac_col <- which(endsWith(colnames(simulations[[1]]), paste0("%%", 0)))
    controls_mat <- sapply(simulations, function(sim){
      sim[,frac_col]
    })
    colnames(controls_mat) <- unname(sapply(simulations, function(sim){
      gsub("%%*.", "", colnames(sim)[frac_col])
    }))


    # Score and transform simulations
    sim_transformed <- scoreTransform(mat = ctoi_mat, signatures, trans_models, is_controls = FALSE)
    controls_mat_uniq <- controls_mat[,!duplicated(colnames(controls_mat))]
    controls_mat_transformed <- scoreTransform(mat = controls_mat_uniq, signatures, trans_models, is_controls = TRUE)
    # Undo unique
    controls_mat_transformed <- sapply(colnames(controls_mat), function(ctrl){
      controls_mat_transformed[,ctrl]
    })
    controls_mat_transformed <- controls_mat_transformed[colnames(sim_transformed), ]

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
#' @param sim_noise description
#' @param ct_sims description
#' @param regGamma description
#' @param sims_sample_frac description
#' @param nCores description
#' @param mix description
#' @param simMethod description
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, data_type, mix = NULL, lineage_file = NULL, weightGenes = TRUE, medianGEP = TRUE, seed = 123, probs = c(0.01, 0.05, 0.1, 0.25, 0.333, 0.49),
                        sim_fracs = c(0, seq(0.01, 0.099, 0.001), seq(0.1, 0.245, 0.005), seq(0.25, 1, 0.05)), diff_vals = c(1, 1.32, 1.585, 2, 3, 4, 5),
                        min_genes = 5, max_genes = 200, return_sigs = FALSE, sigsFile = NULL, minPBcells = 30, minPBsamples = 10,
                        ct_sims = 20, sims_sample_frac = 0.33, simMethod = "ref_thin", sim_noise = NULL, regGamma = 0.8, nCores = 1){


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
    quantiles_matrix <- makeQuantiles(ref_log, labels, probs, dep_list)
    message("Generating signatures...")
    signatures <- createSignatures(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes = weightGenes, ncores = nCores)

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
  simulations <- makeSimulations(ref, labels, mix, gep_mat, cor_mat, dep_list, sim_fracs, sim_method = simMethod, ctoi_samples_frac = sims_sample_frac, n_sims = ct_sims, noise = sim_noise, ncores = nCores, seed2use = seed)
  message("Scoring simulations...")
  simulations_scored <- scoreSimulations(signatures, simulations, nCores)


  # Filter signatures and train RF model
  message("Filtering signatures and training models...")
  models <- trainModels(simulations_scored, regGamma, nCores, seed2use = seed)
  signatures <- signatures[unlist(models$sigs_filtered)]
  models <- models[,-3]


  # Get spillover matrix
  message("Generating spillover matrix...")
  # spill_mat <- getSpillOverMat(simulations, signatures, dep_list, trans_models, n_sims = 1, frac2use = 0.25)
  spill_mat <- matrix()

  # Save results in S4 object
  xCell2Sigs.S4 <- new("xCell2Signatures",
                       signatures = signatures,
                       dependencies = dep_list,
                       models = models,
                       spill_mat = spill_mat,
                       genes_used = rownames(ref))


  return(xCell2Sigs.S4)

}
