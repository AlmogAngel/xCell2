##########################################################################################
# Train signatures
##########################################################################################

# NOTE (!):
#  - No underscore in cell types labels
#  - Cell type onthology with colon in example: CL:0000545 (not CL_0000545)

# USAGE:
# (1) ref - reference expression matrix: genes x samples ("matrix", "dgCMatrix", "Matrix")
# (2) labels - a data frame in which the rows correspond to samples in ref:
#   (a) first column for cell type onthology - named "ont" (character)
#   (b) second column for cell type name - named "label" (character)
#   (c) third column for cell type sample - named "sample" (character)
#   (d) fourth column for cell type dataset - named "dataset" (character)
# (3) data_type - reference data type ("rnaseq", "array", "scrnaseq")
# (4) lineage_file - path to the file generated with xCell2GetLineage and was manually checked (if NULL dependencies list will be made automatically)


library(tidyverse)

# Remove
if (1 == 0) {
  data_type = "sc"; mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02))
  probs = c(.1, .25, .33333333, .5); diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5)
  min_genes = 5; max_genes = 500; is_10x = TRUE
}


xCell2Train <- function(ref, labels, data_type, lineage_file = NULL, mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02)),
                        probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                        min_genes = 5, max_genes = 500, is_10x = FALSE){


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
