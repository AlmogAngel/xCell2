#' xCell2Train function
#'
#' This function generates signatures for each cell type.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import readr
#' @import BiocParallel
#' @importFrom minpack.lm nlsLM
#' @importFrom glmnet cv.glmnet
#' @importFrom Rfast rowMedians rowmeans rowsums Sort
#' @importFrom Matrix rowMeans rowSums colSums
#' @importFrom singscore rankGenes simpleScore
#' @param ref A reference gene expression matrix (genes in rows samples/cells in columns).
#' @param labels A data frame in which the rows correspond to samples in the ref. The data frame must have four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545" or NA if there is no ontology).
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample that match the column name in ref.
#'   "dataset": sample's source dataset or subject (for single-cell).
#' @param ref_type Gene expression data type: "rnaseq" for bulk RNA-Seq, "array" for micro-array, or "sc" for scRNA-Seq.
#' @param seed Set seed for reproducible results (optional).
#' @param min_pb_cells For scRNA-Seq reference only - minimum number of cells in the pseudo-bulk (optional).
#' @param min_pb_samples For scRNA-Seq reference only - minimum number of pseudo-bulk samples (optional).
#' @param min_sc_genes description
#' @param base_count description
#' @param use_ontology A Boolean for using ontological integration (TRUE)
#' @param lineage_file Path to the cell type lineage file generated with `xCell2GetLineage` function (optional).
#' @param probs A numeric vector of probability thresholds to be used for generating signatures (optional).
#' @param diff_vals A numeric vector of delta values to be used for generating signatures (optional).
#' @param top_genes_frac Use for calibration of signatures generation (remove!)
#' @param sim_fracs A vector of mixture fractions to be used in signature filtering (optional).
#' @param min_genes The minimum number of genes to include in the signature (optional).
#' @param max_genes The maximum number of genes to include in the signature (optional).
#' @param sigsFile description
#' @param top_sigs_frac description
#' @param filter_sigs description
#' @param n_sims description
#' @param num_threads description
#' @param human2mouse description
#' @param mix description
#' @param noise_level description
#' @param filtering_data description
#' @param add_essential_genes description
#' @param return_analysis description
#' @param return_signatures description
#' @param use_sillover description
#' @param spillover_alpha description
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref,
                        mix = NULL,
                        labels,
                        ref_type,
                        human2mouse = FALSE,
                        lineage_file = NULL,
                        filtering_data = NULL,
                        seed = 123,
                        num_threads = 1,
                        return_signatures = FALSE,
                        return_analysis = FALSE,
                        use_sillover = TRUE,
                        spillover_alpha = 0.2,
                        # For tuning
                        min_pb_cells = 30,
                        min_pb_samples = 10,
                        min_sc_genes = 1e4,
                        base_count = 3,
                        use_ontology = TRUE,
                        probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4),
                        diff_vals = round(c(log2(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3),
                        min_genes = 3,
                        max_genes = 150,
                        top_genes_frac = 1,
                        filter_sigs = TRUE,
                        sim_fracs = c(0, seq(0.01, 0.25, 0.01)),
                        n_sims = 10,
                        noise_level = 0.2,
                        top_sigs_frac = 0.05,
                        add_essential_genes = TRUE
){

  set.seed(seed)

  # Validate inputs
  inputs_validated <- validateInputs(ref, labels, ref_type)
  ref <- inputs_validated$ref
  labels <- inputs_validated$labels

  # Generate pseudo bulk from scRNA-Seq reference
  if (ref_type == "sc") {
    message("Converting scRNA-seq reference to pseudo bulk...")
    pb_data <- sc2pseudoBulk(ref, labels, min_pb_cells, min_pb_samples)
    ref <- pb_data$ref
    labels <- pb_data$labels
  }

  # Prepare reference and mixture data: human to mouse genes transformation, normalization, base_counts, log2 transformation, shared genes
  out <- prepRefMix(ref, mix, ref_type, min_sc_genes, base_count, human2mouse)
  ref <- out$ref.out
  mix <- out$mix.out
  shared_genes <- rownames(ref)

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  gep_mat <- makeGEPMat(ref, labels)
  cor_mat <- getCellTypeCorrelation(gep_mat, ref_type)

  # Get cell type dependencies list
  if (use_ontology) {
    message("Loading dependencies...")
    if (is.null(lineage_file)) {
      dep_list <- xCell2::xCell2GetLineage(labels, out_file = NULL)
    }else{
      dep_list <- getDependencies(lineage_file)
    }
  }else{
    message("Skipping ontological integration")
    dep_list <- NULL
  }

  # Generate signatures
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs, num_threads)
  message("Generating signatures...")
  signatures <- createSignatures(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, top_genes_frac, num_threads)

  if (return_signatures) {
    xCell2.S4 <- new("xCell2Object",
                     signatures = signatures,
                     dependencies = dep_list,
                     params = NULL,
                     spill_mat = NULL,
                     genes_used = shared_genes)
    return(xCell2.S4)
  }

  if (filter_sigs) {
    # Generate simulations
    message("Generating simulations...")
    simulations <- makeSimulations(ref, mix, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims, noise_level, num_threads)

    # Filter signatures
    message("Filtering signatures...")
    signatures <- filterSignatures(shared_genes, labels, filtering_data, simulations, signatures, top_sigs_frac, add_essential_genes, human2mouse, num_threads)
  }

  # Learn linear transformation parameters
  message("Learning linear transformation and spillover parameters...")
  params <- learnParams(gep_mat, cor_mat, signatures, dep_list, ref_type, sim_fracs, frac2use = 0.25, num_threads)
  spill_mat <- params$spillmat
  params <- params$params


  # Save results in S4 object
  xCell2.S4 <- new("xCell2Object",
                   signatures = signatures,
                   dependencies = dep_list,
                   params = params,
                   spill_mat = spill_mat,
                   genes_used = shared_genes)

  message("Your custom xCell2 reference object is ready!")
  message("> Please consider sharing your xCell2 reference with others here: https://dviraran.github.io/xCell2ref")


  if (return_analysis) {
    message("Running xCell2Analysis...")
    res <- xCell2::xCell2Analysis(mix, xcell2object = xCell2.S4, spillover = use_sillover, spillover_alpha = spillover_alpha, num_threads)
    return(res)
  }else{
    return(xCell2.S4)
  }

}
