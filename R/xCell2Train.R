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
prepareRef <- function(ref, data_type, bulk_pseudo_count){

  if (data_type == "sc") {
    if(max(ref) >= 50){
      message("Normalizing and transforming scRNA-Seq reference to log1p-space (maximum expression value >= 50).")
      genes_names <- rownames(ref)
      ref.srt <- Seurat::CreateSeuratObject(counts = ref)
      ref.srt <- Seurat::NormalizeData(ref.srt, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
      ref.norm <-  ref.srt@assays$RNA@data

      # Check if Seurat changed genes names
      if (!all(rownames(ref.norm) %in% genes_names)) {
        stop("Seurat genes name error")
      }


      return(ref.norm)
    }else{
      message("Assuming reference is already in log1p-space (maximum expression value < 50).")
      return(ref)
    }

  }else{

    if(max(ref) >= 50){
      message("Transforming reference to log2-space (maximum expression value >= 50).")
      ref.norm <- log2(ref+bulk_pseudo_count)
      return(ref.norm)
    }else{
      message("Assuming reference is already in log2-space (maximum expression value < 50).")
      return(ref)
    }

  }

}
sc2pseudoBulk <- function(ref, labels, min_n_cells, min_groups, seed = 123){

  set.seed(seed)

  celltypes <- unique(labels$label)

  groups_list <- lapply(celltypes, function(ctoi){

    ctoi_samples <- labels[labels$label == ctoi,]$sample

    # Calculate maximum possible number of groups given min_n_cells
    num_groups <- ceiling(length(ctoi_samples) / min_n_cells)
    if (num_groups < min_groups) {
      num_groups <- min_groups
    }

    # Generate min_groups pseudo samples of CTOI
    if (length(ctoi_samples) > min_groups) {

      ctoi_samples_shuffled <- sample(ctoi_samples, length(ctoi_samples))
      list_of_ctoi_samples_shuffled <- split(ctoi_samples_shuffled, ceiling(seq_along(ctoi_samples_shuffled) / (length(ctoi_samples_shuffled) / num_groups)))

      sapply(list_of_ctoi_samples_shuffled, function(ctoi_group){
        if (length(ctoi_group) == 1) {
          ref[,ctoi_group]
        }else{
          if("matrix" %in% class(ref)) Rfast::rowmeans(ref[,ctoi_group]) else Matrix::rowMeans(ref[,ctoi_group])
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
makePureCTMat <- function(ref, labels, use_median){

  celltypes <- unique(labels$label)

  pure_ct_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      if(use_median){
        type_vec <- if("matrix" %in% class(ref)) Rfast::rowMedians(ref[,type_samples]) else sparseMatrixStats::rowMedians(ref[,type_samples])
      }else{
        type_vec <- if("matrix" %in% class(ref)) Rfast::rowmeans(ref[,type_samples]) else Matrix::rowMeans(ref[,type_samples])
      }
    }
  })
  rownames(pure_ct_mat) <- rownames(ref)

  return(pure_ct_mat)
}
getCellTypeCorrelation <- function(pure_ct_mat, data_type){

  celltypes <- colnames(pure_ct_mat)

  if (data_type != "sc") {

    # Use top 10% most variable genes
    genes_var <- apply(pure_ct_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.9, na.rm=TRUE)
    pure_ct_mat <- pure_ct_mat[genes_var > most_var_genes_cutoff,]

  }else{

    # Use top 1% most variable genes
    genes_var <- apply(pure_ct_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.99, na.rm=TRUE)
    pure_ct_mat <- pure_ct_mat[genes_var > most_var_genes_cutoff,]

  }


  # Make correlation matrix
  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

  # TODO: Change for loop to apply function to measure time
  for (i in 1:nrow(lower_tri_coord)) {
    celltype_i <- rownames(cor_mat)[lower_tri_coord[i, 1]]
    celltype_j <- colnames(cor_mat)[lower_tri_coord[i, 2]]
    cor_mat[lower_tri_coord[i, 1], lower_tri_coord[i, 2]] <- cor(pure_ct_mat[,celltype_i], pure_ct_mat[,celltype_j], method = "spearman")
    cor_mat[lower_tri_coord[i, 2], lower_tri_coord[i, 1]] <- cor(pure_ct_mat[,celltype_i], pure_ct_mat[,celltype_j], method = "spearman")
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
makeQuantiles <- function(ref, labels, probs, dep_list, include_descendants){

  celltypes <- unique(labels[,2])

  quantiles_mat_list <-  pbapply::pblapply(celltypes, function(type){

    if (include_descendants) {
      # Include all the descendants of the cell type for quantiles calculations
      descen_cells <- dep_list[[type]]$descendants
      type_samples <- labels[,2] == type | labels[,2] %in% descen_cells
    }else{
      type_samples <- labels[,2] == type
    }


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
# TODO: Check weight_genes = FALSE
createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Set signature parameters
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)

    # Generate signatures
    type_sigs <- list()
    for(i in 1:nrow(param.df)){

      # Get a Boolean matrices with genes that pass the quantiles criteria
      diff <- param.df[i, ]$diff_vals # difference criteria
      lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile criteria
      upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile criteria
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


      # Round and sort genes scores
      top_scores <- sort(unique(round(gene_passed-0.5)), decreasing = TRUE)

      # Take top 10 highest scores
      # TODO: test different n_tops parameters (currently 10 is default)
      n_top <- ifelse(length(top_scores) > 10, 10, length(top_scores))

      for (score in top_scores[1:n_top]) {

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

    if (length(type_sigs) == 0) {
      warning(paste0("No signatures found for ", type))
    }


    # Remove duplicate signatures
    type_sigs_sorted <- lapply(type_sigs, function(x) sort(x))
    type_sigs_sorted_collapsed <- sapply(type_sigs_sorted, paste, collapse = ",")
    duplicated_sigs <- duplicated(type_sigs_sorted_collapsed)
    type_sigs <- type_sigs[!duplicated_sigs]


    return(type_sigs)
  }


  celltypes <- unique(labels[,2])

  all_sigs <- pbapply::pblapply(celltypes, function(type){
    getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes)
  })
  all_sigs <- unlist(all_sigs, recursive = FALSE)


  if (length(all_sigs) == 0) {
    warning("No signatures found for reference!")
  }


  return(all_sigs)
}
makeSimulations <- function(ref, labels, pure_ct_mat, cor_mat, dep_list, sim_fracs, ctoi_samples_frac = 1, n_sims = 1){

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
  makeFractionMatrix <- function(ref, sim_fracs, samples2use, control){

    if (length(samples2use) == 1) {
      mean_expression <- ref[,samples2use]
    }else{
      mean_expression <- if("matrix" %in% class(ref)) Rfast::rowmeans(ref[,samples2use]) else Matrix::rowMeans(ref[,samples2use])
    }

    if (control) {
      sim_fracs <- 1-sim_fracs
    }

    frac_mat <- matrix(rep(mean_expression, length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs)) %*% diag(sim_fracs)
    rownames(frac_mat) <- rownames(ref)

    return(frac_mat)
  }

  sim_list <- pbapply::pblapply(celltypes, function(ctoi){

    # Generate CTOI samples pool
    ctoi_samples_pool <- makeSamplesPool(labels, ctoi)

    # Number of CTOI samples to use
    n_samples_sim <- round(length(ctoi_samples_pool) * ctoi_samples_frac)
    n_samples_sim <- ifelse(n_samples_sim < 1, 1, n_samples_sim)


    # Get control cell types
    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
    controls <- celltypes[!celltypes %in% dep_cts]
    control <- names(sort(cor_mat[ctoi, controls])[1])


    # Generate n_sims simulations
    ctoi_sim_list <- lapply(1:n_sims, function(i){

      # Make CTOI fraction matrix
      samples2use <- ctoi_samples_pool[1:n_samples_sim]
      ctoi_frac_mat <- makeFractionMatrix(ref, sim_fracs, samples2use, control = FALSE)

      # Move samples2use to the end of the vector
      ctoi_samples_pool <- c(ctoi_samples_pool[!ctoi_samples_pool %in% samples2use], samples2use)

      # Make control(s) fractions matrix
      samples2use <- labels[labels$label == control,]$sample
      control_frac_mat <- makeFractionMatrix(ref, sim_fracs, samples2use, control = TRUE)


      # Combine CTOI and control(s) fractions matrix
      simulation <- ctoi_frac_mat + control_frac_mat
      colnames(simulation) <- paste0(control, "%%", sim_fracs)

      simulation
    })

    if (n_sims == 1) {
      ctoi_sim_list[[1]]
    }else{
      names(ctoi_sim_list) <- paste0("sim-", 1:n_sims)
      ctoi_sim_list
    }
  })
  names(sim_list) <- celltypes

  return(sim_list)

}
filterSignatures <- function(simulations_scored){


  top_cor_sigs <- bind_rows(simulations_scored) %>%
    mutate(sim_frac = as.numeric(sim_frac)) %>%
    group_by(celltype, signature, sim_id) %>%
    summarise(cor = cor(sim_frac, score)) %>%
    summarise(cor = mean(cor)) %>%
    top_n(n = max(50, 0.5*n()), wt = cor) %>%
    pull(signature)

  return(top_cor_sigs)


}
scoreSimulations <- function(signatures, simulations, dep_list, n_sims){


  celltypes <- names(simulations)
  sims_scored <- pbapply::pblapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]

    if (n_sims == 1) {
      sim_ranked <- singscore::rankGenes(simulations[[ctoi]])

      scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      }))
      colnames(scores) <- colnames(sim_ranked)

      as_tibble(scores, rownames = "signature") %>%
        pivot_longer(cols = -signature, values_to = "score") %>%
        separate(name, into = c("control", "frac"), sep = "%%") %>%
        separate(signature, into = "celltype", sep = "#", remove = FALSE, extra = "drop")

    }else{
      lapply(simulations[[ctoi]], function(sim){

        sim_ranked <- singscore::rankGenes(sim)

        scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
          singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
        }))
        colnames(scores) <- colnames(sim)

        as_tibble(scores, rownames = "signature") %>%
          pivot_longer(cols = -signature, values_to = "score") %>%
          separate(name, into = c("control", "frac"), sep = "%%") %>%
          separate(signature, into = "celltype", sep = "#", remove = FALSE, extra = "drop")

      }) %>%
        bind_rows(., .id = "sim_id")
    }

  })
  names(sims_scored) <- celltypes

  return(sims_scored)

}
trainModels <- function(simulations_scored, params, modelType, seed){

  set.seed(seed)

  fitModel <- function(data, modelParams = params, model_type = modelType){


    train_mat <- data %>%
      select(signature, frac, score) %>%
      mutate(frac = as.numeric(frac)) %>%
      pivot_wider(names_from = signature, values_from = score) %>%
      as.matrix()

    if (model_type == "rf") {


      gamma <- modelParams$gamma
      mtry2use <- ncol(train_mat[,-1]) / modelParams$mtry

      RF <- RRF::RRF(train_mat[,-1], train_mat[,1], flagReg = 0, importance = TRUE, mtry = mtry2use, ntree = modelParams$ntree, nodesize = modelParams$nodesize)
      RF_imp <- RF$importance[,"%IncMSE"] / max(RF$importance[,"%IncMSE"])
      RRF <- RRF::RRF(train_mat[,-1], train_mat[,1], flagReg = 1, coefReg = (1-gamma) + gamma*RF_imp)
      return(RRF)
    }


    if (model_type == "xgb") {


      if (XGBparams$objective == "reg:gamma") {
        train_mat[,1][train_mat[,1] == 0] <- 0.000000001
      }

      train_mat <- xgboost::xgb.DMatrix(data = train_mat[,-1], label = train_mat[,1])

      model <- xgboost::xgb.train(
        params = modelParams,
        data = train_mat,
        nrounds = 100
      )

      return(model)
    }

  }

  enframe(simulations_scored, name = "celltype", value = "data") %>%
    rowwise() %>%
    mutate(model = list(fitModel(data))) %>%
    dplyr::select(celltype, model) %>%
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
#' @slot transformation_models data frame of cell type transformation models
#' @slot spill_mat matrix of cell types spillover
#' @slot genes_used character vector of genes names used to train the signatures
#' @importFrom methods new
# Create S4 object for the new reference
setClass("xCell2Signatures", slots = list(
  signatures = "list",
  dependencies = "list",
  transformation_models = "data.frame",
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
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures VariableFeatures
#' @importFrom Rfast rowMedians rowmeans
#' @importFrom pbapply pblapply pbsapply
#' @importFrom sparseMatrixStats rowMedians
#' @importFrom RRF RRF
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom Matrix rowMeans
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
#' @param params description
#' @param sigsFile description
#' @param filter_sigs description
#' @param modelType description
#' @param bulkPseudoCount description
#' @param minPBcells description
#' @param minPBgroups description
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, data_type, lineage_file = NULL, weightGenes = FALSE,
                        sim_fracs = c(0, 0.001, 0.002, 0.004, 0.006, 0.008, seq(0.01, 1, 0.01)), diff_vals = c(1, 1.32, 1.585, 2, 3, 4, 5),
                        min_genes = 5, max_genes = 200, filter_sigs = TRUE, sigsFile = NULL, params = list(), modelType = "rf", bulkPseudoCount = 1, minPBcells = 30, minPBgroups = 10){


  # Validate inputs
  inputs_validated <- validateInputs(ref, labels, data_type)
  ref <- inputs_validated$ref
  labels <- inputs_validated$labels

  # Prepare reference
  ref <- prepareRef(ref, data_type, bulk_pseudo_count = bulkPseudoCount)

  # Generate pseudo bulk from scRNA-Seq reference
  if (data_type == "sc") {
    message("Converting scRNA-seq reference to pseudo bulk...")
    ps_data <- sc2pseudoBulk(ref, labels, min_n_cells = minPBcells, min_groups = minPBgroups)
    ref <- ps_data$ref
    labels <- ps_data$labels
  }

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  pure_ct_mat <- makePureCTMat(ref, labels, use_median = TRUE)
  cor_mat <- getCellTypeCorrelation(pure_ct_mat, data_type)

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
    if (data_type != "sc") {
      probs <- c(0.01, 0.05, 0.1, 0.25, 0.333, 0.49)
    }else{
      probs <- c(0.1, 0.15, 0.2, 0.25, 0.333, 0.49)
      # Adjust diff values to log1p
      fold_change_vals <- round(2^diff_vals, 4)
      diff_vals <- round(log1p(fold_change_vals - 1), 3)
    }
    message("Calculating quantiles...")
    quantiles_matrix <- makeQuantiles(ref, labels, probs, dep_list, include_descendants = FALSE)
    message("Generating signatures...")
    signatures <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes = weightGenes)

    if (!filter_sigs) {
      return(signatures)
    }
  }else{
    message("Loading signatures...")
    signatures <- readRDS(sigsFile)
  }


  # Make simulations
  message("Generating simulations...")
  simulations <- makeSimulations(ref, labels, pure_ct_mat, cor_mat, dep_list, sim_fracs, ctoi_samples_frac = 1, n_sims = 1)
  simulations_scored <- scoreSimulations(signatures, simulations, dep_list, n_sims = 1)


  # TODO: Filter signatures
  message("Filtering signatures...")


  # Get transformation models
  message("Training models...")
  trans_models <- trainModels(simulations_scored, params, modelType, seed = 123)


  # Get spillover matrix
  message("Generating spillover matrix...")
  spill_mat <- getSpillOverMat(simulations, signatures, dep_list, trans_models, n_sims = 1, frac2use = 0.25)

  # Save results in S4 object
  xCell2Sigs.S4 <- new("xCell2Signatures",
                       signatures = signatures,
                       dependencies = dep_list,
                       transformation_models = trans_models,
                       spill_mat = spill_mat,
                       genes_used = rownames(ref))

  return(xCell2Sigs.S4)

}
