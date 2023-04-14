#' xCell2Train function
#'
#' This function generates signatures for each cell type.
#'
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


  source("utils.R")

  # Use only most variable genes for single-cell data
  if (data_type == "sc") {
    message("Looking for most variable genes with Seurat...")
    topVarGenes <- getTopVariableGenes(ref, min_genes = 10000, sensitivity = 15)
    ref <- ref[rownames(ref) %in% topVarGenes,]
  }

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  pure_ct_mat <- makePureCTMat(ref, labels)
  cor_mat <- getCellTypeCorrelation(pure_ct_mat, data_type)

  # Get cell type dependencies list
  message("Loading dependencies...")
  if (is.null(lineage_file)) {
    dep_list <- xCell2GetLineage(labels = labels[,1:2], out_file = NULL)
  }else{
    dep_list <- getDependencies(lineage_file)
  }

  # Generate signatures for each cell type
  message("Generating signatures...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs)
  signatures_collection <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes)

  # Filter signatures
  message("Filtering signatures...")
  filter_signature_out <- filterSignatures(ref, labels, pure_ct_mat, dep_list, signatures_collection, mixture_fractions, grubbs_cutoff = 0.8, simulations_cutoff = 0.8)
  scores_mat_pure_tidy <- filter_signature_out$scoreMatTidy
  signatures_collection_filtered <- filter_signature_out$sigCollectionFilt
  # plotHeatMap("Neutrophils", scores_mat_pure_tidy, signatures_collection_filtered = NULL, cor_mat)
  # plotHeatMap("Neutrophils", scores_mat_pure_tidy, signatures_collection_filtered = signatures_collection_filtered, cor_mat)

  # TODO: Weight signatures with Elastic Net
  # source("train_models.R")
  # source("R/train_models_tmp.R")
  # models <- trainModels(ref, labels, dep_list, pure_ct_mat_test, signatures_collection_filtered, mixture_fractions)

  # TODO: Linear tranformation


  # Create S4 object for the new reference
  setClass("xCell2 Signatures", slots = list(
    labels = "data.frame",
    dependencies = "list",
    all_signatures = "GeneSetCollection",
    filtered_signatures = "GeneSetCollection"
  ))


  xCell2Ref.S4 <- new("xCell2 Signatures", labels = labels, dependencies = dep_list, all_signatures = signatures_collection, filtered_signatures = signatures_collection_filtered)

  return(xCell2Ref.S4)

}
