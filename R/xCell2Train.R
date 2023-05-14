validateInputs <- function(ref, labels, data_type){
  if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
    stop("ref should be as matrix.")
  }

  if (!"data.frame" %in% class(labels)) {
    stop("labels should be as dataframe.")
  }

  if (!data_type %in% c("rnaseq", "array", "sc")) {
    stop("data_type should be rnaseq, array or scrnaseq.")
  }

  if (sum(grepl("_", labels$label)) != 0 | sum(grepl("_", rownames(ref))) != 0) {
    message("Changing underscores to dashes in genes / cell-types labels!")
    labels$label <- gsub("_", "-", labels$label)
    rownames(ref) <- gsub("_", "-", rownames(ref))
  }

}
getTopVariableGenes <- function(ref, min_genes){

  ref.srt <- Seurat::CreateSeuratObject(counts = ref)
  ref.srt <- Seurat::FindVariableFeatures(ref.srt, selection.method = "vst")
  plot1 <- Seurat::VariableFeaturePlot(ref.srt)
  genesVar <- plot1$data$variance.standardized
  names(genesVar) <- rownames(plot1$data)
  genesVar <- sort(genesVar, decreasing = TRUE)

  if (length(genesVar) < min_genes) {
    genesVar <- names(genesVar[1:min_genes])
  }

  return(genesVar)
}
makePureCTMat <- function(ref, labels){

  celltypes <- unique(labels$label)

  pure_ct_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      type_vec <- if("matrix" %in% class(ref)) Rfast::rowMedians(ref[,type_samples]) else sparseMatrixStats::rowMedians(ref[,type_samples])
    }
  })
  rownames(pure_ct_mat) <- rownames(ref)

  return(pure_ct_mat)
}
getCellTypeCorrelation <- function(pure_ct_mat, data_type){

  celltypes <- colnames(pure_ct_mat)

  # Use 50% most expressed genes because single-cell is zero inflated
  if (data_type != "sc") {
    mean_gene_expression <- Rfast::rowmeans(pure_ct_mat)
    high_gene_expression_cutoff <- quantile(mean_gene_expression, 0.5, na.rm=TRUE)
    top_expressed_gene <- mean_gene_expression > high_gene_expression_cutoff
    pure_ct_mat <- pure_ct_mat[top_expressed_gene,]
  }

  # Use top 85% most variable genes
  genes_var <- apply(pure_ct_mat, 1, var)
  most_var_genes_cutoff <- quantile(genes_var, 0.85, na.rm=TRUE)
  pure_ct_mat <- pure_ct_mat[genes_var > most_var_genes_cutoff,]

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

  quantiles_matrix <-  pbapply::pblapply(celltypes, function(type){

    # Include all the descendants of the cell type in the quantiles calculations
    if (include_descendants) {
      descen_cells <- dep_list[[type]]$descendants
      type_samples <- labels[,2] == type | labels[,2] %in% descen_cells
    }else{
      type_samples <- labels[,2] == type
    }

    # If there is one sample for this cell type -> duplicate the sample to make a data frame
    if (sum(type_samples) == 1) {
      type.df <- cbind(ref[,type_samples], ref[,type_samples])
    }else{
      type.df <- ref[,type_samples]
    }

    # Calculate quantiles
    quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1-probs))), na.rm=TRUE))
  })
  names(quantiles_matrix) <- celltypes

  return(quantiles_matrix)
}
createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Weight cell types by correlation
    type_weights <- cor_mat[type, not_dep_celltypes]

    # Signature parameters
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)

    # Generate signatures
    type_sigs <- list()
    for(i in 1:nrow(param.df)){

      # Get a Boolean matrices with genes that pass the quantiles criteria
      diff <- param.df[i, ]$diff_vals # diffrence criteria
      lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile criteria
      upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile criteria

      diff_genes.mat <- sapply(not_dep_celltypes, function(x){
        get(type, quantiles_matrix)[lower_prob,] > get(x, quantiles_matrix)[upper_prob,] + diff
      })

      # Score genes using weights
      gene_scores <- apply(diff_genes.mat, 1, function(x){
        sum(type_weights[which(x)])
      })

      gene_passed <- gene_scores[gene_scores > 0]

      # If less than min_genes passed move to next parameters
      if (length(gene_passed) < min_genes) {
        next
      }

      # Save signatures
      gene_passed <- sort(gene_passed, decreasing = TRUE)
      for (n_genes in round(seq(from = min_genes, to = max_genes, length.out = 8))) {

        if (length(gene_passed) < n_genes) {
          break
        }

        sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, sep = "_")
        type_sigs[[sig_name]] <- GSEABase::GeneSet(names(gene_passed[1:n_genes]), setName = sig_name)
      }

    }

    if (length(type_sigs) == 0) {
      warning(paste0("No signatures found for ", type))
    }

    return(type_sigs)
  }

  celltypes <- unique(labels[,2])

  all_sigs <- pbapply::pblapply(celltypes, function(type){
    getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes)
  })
  all_sigs <- unlist(all_sigs)


  if (length(all_sigs) == 0) {
    warning("No signatures found for reference!")
  }

  # Make GeneSetCollection object
  signatures_collection <- GSEABase::GeneSetCollection(all_sigs)


  return(signatures_collection)
}
filterSignatures <- function(ref, labels, mixture_fractions, dep_list, signatures_collection, top_sigs = 0.05){

  # Make simulations
  makeSimulations <- function(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = TRUE, seed = 123){

    set.seed(seed)

    makeFractionMatrixCTOI <- function(ref, mixture_fractions, ctoi_samples_pool){


      # Make CTOI fraction matrix
      ctoi_frac_mat <- matrix(NA, nrow = nrow(ref), ncol = length(mixture_fractions), dimnames = list(rownames(ref), mixture_fractions))
      for (i in 1:length(mixture_fractions)) {
        frac <- mixture_fractions[i]
        ctoi_samples2use <- ctoi_samples_pool[1:length(mixture_fractions)] # Choose the first samples (homogeneous by datasets)
        frac_fracs <- diff(c(0, sort(runif(length(ctoi_samples2use)-1, min = 0, max = frac)), frac)) # Generate random fraction for each sample to sum to frac in mixture_fractions
        ctoi_frac_mat[,i] <- Rfast::rowsums(ref[,ctoi_samples2use] %*% diag(frac_fracs)) # Sum all fractions
        ctoi_samples_pool <- c(ctoi_samples_pool[!ctoi_samples_pool %in% ctoi_samples2use], ctoi_samples2use) # Move ctoi_samples2use to be last
      }

      out <- list(frac_mat = ctoi_frac_mat,
                  samples_pool = ctoi_samples_pool)

      return(out)
    }

    makeFractionMatrixControls <- function(ref, labels, ctoi, dep_cts, mixture_fractions){


      # Pick controls
      control_samples <- labels %>%
        filter(!label %in% dep_cts) %>%
        slice_sample(n=1, by = label) %>%
        slice_sample(n=length(mixture_fractions)) %>%
        pull(sample)
      controls_expression <- ref[,control_samples]

      controls_mat <- sapply(1-mixture_fractions, function(frac){
        random_fracs <- diff(c(0, sort(runif(ncol(controls_expression)-1, min = 0, max = frac)), frac)) # Generate random numbers from the controls from a uniform distribution that sum to frac
        Rfast::rowsums(controls_expression %*% diag(random_fracs))
      })

      return(controls_mat)

    }


    celltypes <- unique(labels[,2])
    ctoi_sim_list <- list()
    sim_list <- pbapply::pblapply(celltypes, function(ctoi){

      # Sort CTOI samples to be homogeneous by datasets
      ctoi_samples_pool <- c()
      while(!all(labels[labels$label == ctoi,]$sample %in% ctoi_samples_pool)) {
        ctoi_samples_pool <- c(ctoi_samples_pool,
                               labels %>%
                                 filter(label == ctoi & !sample %in% ctoi_samples_pool) %>%
                                 slice_head(n = 1, by = dataset) %>%
                                 pull(sample))
      }

      if (length(ctoi_samples_pool) < length(mixture_fractions)) {
        ctoi_samples_pool <- rep(ctoi_samples_pool, length(mixture_fractions))
      }

      dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
      n_ctoi_ds <- length(unique(labels[labels$label == ctoi,]$dataset))


      for (i in 1:n_ct_sim) {

        ctoi_frac_mat_results <- makeFractionMatrixCTOI(ref, mixture_fractions, ctoi_samples_pool)
        ctoi_frac_mat <- ctoi_frac_mat_results$frac_mat
        ctoi_samples_pool <- ctoi_frac_mat_results$samples_pool

        # Make Controls fraction matrix
        controls_frac_mat <- makeFractionMatrixControls(ref, labels, ctoi, dep_cts, mixture_fractions)

        # Combine CTOI and controls fractions matrix
        simulation <- ctoi_frac_mat + controls_frac_mat

        # Add noise
        if (add_noise) {
          noise_sd <- 1 + 1/n_ctoi_ds
          noise <- matrix(rnorm(nrow(simulation) * ncol(simulation), mean = 0, sd = noise_sd),
                          nrow = nrow(simulation), ncol = ncol(simulation))
          simulation <- simulation + noise
          simulation <- pmax(simulation, 0)
        }

        ctoi_sim_list[[paste0("sim-", i)]] <- simulation
      }

      ctoi_sim_list
    })
    names(sim_list) <- celltypes

    return(sim_list)

  }
  sim_list <- makeSimulations(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = TRUE, seed = 123)

  # Score simulations
  scoreSimulations <- function(signatures_collection, sim_list){

    celltypes <- names(sim_list)
    out <- pbapply::pblapply(celltypes, function(ctoi){
      signatures_ctoi <- signatures_collection[startsWith(names(signatures_collection), paste0(ctoi, "#"))]

      scores_ctoi <- lapply(sim_list[[ctoi]], function(sim){
        sim_ranked <- singscore::rankGenes(sim)
        scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
          singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
        }))

        rownames(scores) <- names(signatures_ctoi)
        colnames(scores) <- colnames(sim)

        scores
      })
      scores_ctoi

    })
    names(out) <- celltypes
    return(out)
  }
  scores_list <- scoreSimulations(signatures_collection, sim_list)

  # Get simulation correlations results
  getSimulationCor <- function(scores){
    fractions <- as.numeric(colnames(scores))
    sig_cors <- apply(scores, 1, function(sig){
      cor(sig, fractions, method = "spearman")
    })
    return(sig_cors)
  }
  sim_scores_cors_tidy <- enframe(scores_list, name = "celltype") %>%
    unnest_longer(value, indices_to = "sim_id", values_to = "scores") %>%
    rowwise() %>%
    mutate(cor = list(getSimulationCor(scores))) %>%
    unnest_longer(cor, indices_to = "signature") %>%
    ungroup()

  # Filter by top signatures
  top_sigs <- sim_scores_cors_tidy %>%
    group_by(celltype, signature) %>%
    summarise(mean_cor = mean(cor)) %>%
    group_by(celltype) %>%
    top_frac(n=top_sigs, wt=mean_cor) %>%
    pull(signature)
  signatures_filtered <- signatures_collection[names(signatures_collection) %in% top_sigs]

  # Filter by Grubbs' test (?)

  # Return: (1) filtered_signatures (2) simulations + signatures scores
  out <- list(signatures_filtered = signatures_filtered,
              sim_scores_cors = sim_scores_cors_tidy)

  return(out)
}
# TODO:
transformScores <- function(scores_mat_tidy, signatures_filtered){

  # TODO: check mixture_fraction hard numbers in the filters

  # Find calibration values
  calibration_values <- scores_mat_tidy %>%
    filter(signature %in% pull(signatures_filtered, signature) & signature_ct == mixture_ct & mixture_fraction < 1) %>%
    # Take the mean of all signatures' scores per cell type per mixture fraction
    group_by(signature_ct, mixture_fraction) %>%
    summarise(mean_score = mean(score)) %>%
    # The shift value is the minimum score per cell type
    mutate(shift_value = min(mean_score)) %>%
    # The scale factor is the slope of the new linear line (i.e., 0.25/shifted_score(25%))
    mutate(shifted_score = mean_score - shift_value) %>%
    mutate(scale_factor = max(as.numeric(mixture_fraction))/max(shifted_score)) %>%
    select(signature_ct, shift_value, scale_factor) %>%
    unique()

  # Transform scores
  scores_transformed <- scores_mat_tidy %>%
    filter(signature %in% pull(signatures_filtered, signature) & mixture_fraction == 0.25) %>%
    group_by(signature_ct, mixture_ct) %>%
    summarise(mean_score = mean(score)) %>%
    left_join(calibration_values, by = "signature_ct") %>%
    rowwise() %>%
    mutate(transformed_score = (mean_score-shift_value)*scale_factor)

  return(scores_transformed)

}


#' @slot labels ...
#' @slot dependencies ...
#' @slot all_signatures ...
#' @slot filtered_signatures ...
#' @name xCell2Signatures
#' @importFrom methods new
# Create S4 object for the new reference
setClass("xCell2Signatures", slots = list(
  labels = "data.frame",
  dependencies = "list",
  all_signatures = "GeneSetCollection",
  filtered_signatures = "GeneSetCollection"
))

# Remove - for debugging
if (0 == 1) {
  data_type = "rnaseq"; lineage_file = NULL; mixture_fractions = c(seq(0, 0.25, 0.03), 1)
  probs = c(.1, .25, .33333333, .5); diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5); min_genes = 5; max_genes = 200
}


#' xCell2Train function
#'
#' This function generates signatures for each cell type.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import Seurat
#' @import Rfast
#' @import pbapply
#' @import sparseMatrixStats
#' @importFrom  GSEABase GeneSetCollection
#' @importFrom  GSEABase GeneSet
#' @import singscore
#' @param ref A reference gene expression matrix.
#' @param labels A data frame in which the rows correspond to samples in the ref. The data frame must have four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545").
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample should match the column name in ref.
#'   "dataset": the cell type sample dataset or subject (for single-cell) as a character.
#' @param data_type Gene expression data type: "rnaseq", "array", or "sc".
#' @param lineage_file (Optional) Path to the cell type lineage file generated with `xCell2GetLineage` function.
#' @param mixture_fractions A vector of mixture fractions to be used in signature filtering.
#' @param probs A vector of probability thresholds to be used for generating signatures.
#' @param diff_vals A vector of delta values to be used for generating signatures.
#' @param min_genes The minimum number of genes to include in the signature.
#' @param max_genes The maximum number of genes to include in the signature.
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, data_type, lineage_file = NULL, mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02)),
                        probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                        min_genes = 5, max_genes = 500){


  # Validate inputs
  validateInputs(ref, labels, data_type)

  # Use only most variable genes for single-cell data
  if (data_type == "sc") {
    message("Looking for most variable genes with Seurat...")
    topVarGenes <- getTopVariableGenes(ref, min_genes = 10000)
    ref <- ref[rownames(ref) %in% topVarGenes,]
  }

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  pure_ct_mat <- makePureCTMat(ref, labels)
  cor_mat <- getCellTypeCorrelation(pure_ct_mat, data_type)

  # Get cell type dependencies list
  message("Loading dependencies...")
  if (is.null(lineage_file)) {
    dep_list <- xCell2::xCell2GetLineage(labels = labels[,1:2], out_file = NULL)
  }else{
    dep_list <- getDependencies(lineage_file)
  }

  # Generate signatures for each cell type
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs, dep_list, include_descendants = TRUE)
  message("Generating signatures...")
  signatures_collection <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes)

  # Filter signatures
  message("Filtering signatures...")
  #filter_signature_out <- filterSignatures(ref, labels, pure_ct_mat, dep_list, signatures_collection, mixture_fractions, grubbs_cutoff = 0.8, simulations_cutoff = 0.8)
  #scores_mat_pure_tidy <- filter_signature_out$scoreMatTidy
  #signatures_collection_filtered <- filter_signature_out$sigCollectionFilt


  # TODO: Weight signatures


  # TODO: Linear transformation


  # TODO: Spillover correction


  xCell2Ref.S4 <- new("xCell2Signatures", labels = labels, dependencies = dep_list, all_signatures = signatures_collection, filtered_signatures = signatures_collection_filtered)

  return(xCell2Ref.S4)

}
