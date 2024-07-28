#' xCell2Analysis function
#'
#' This function takes a matrix of mix gene expression data and a `xCell2Signatures` object containing a set of signatures as input. It performs downstream analysis to identify enriched cell types in the mix sample.
#'
#' @importFrom singscore rankGenes simpleScore
#' @importFrom randomForestSRC predict.rfsrc
#' @import dplyr
#' @import tibble
#' @import purrr
#' @importFrom  pracma lsqlincon
#' @importFrom xgboost xgb.DMatrix
#' @param mix a matrix containing gene expression data
#' @param xcell2object S4 object of `xCell2Object`
#' @param min_shared_genes description
#' @param ref_is_sc description
#' @param raw_scores description
#' @param spillover Boolean - should we use spillover correction on the transformed scores?
#' @param spillover_alpha description
#' @param num_threads description

#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
#' @export
xCell2Analysis <- function(mix, xcell2object, ref_is_sc, min_shared_genes = 0.9, raw_scores = FALSE, spillover = TRUE, spillover_alpha = 0.2, num_threads = 1){

  scoreMix <- function(ctoi, mix_ranked, signatures_ctoi){

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
    })

    rownames(scores) <- colnames(mix_ranked)

    return(scores)
  }


  param <- BiocParallel::MulticoreParam(workers = num_threads)


  # Check reference/mixture genes intersection
  genes_intersect_frac <- length(intersect(rownames(mix), xcell2object@genes_used)) / length(xcell2object@genes_used)
  if (genes_intersect_frac < min_shared_genes) {
    stop("This xCell2 reference shares ", genes_intersect_frac, " genes with the mixtures and min_shared_genes = ", min_shared_genes, ".",
         "\n", "Consider training a new xCell2 reference or adjusting min_shared_genes.")
  }


  # Rank mix gene expression matrix
  mix_ranked <- singscore::rankGenes(mix[xcell2object@genes_used,])

  # Score and predict
  sigs_celltypes <- unique(unlist(lapply(names(xcell2object@signatures), function(x){strsplit(x, "#")[[1]][1]})))

  # Get raw enrichment scores
  res_raw <- BiocParallel::bplapply(sigs_celltypes, function(ctoi){

    signatures_ctoi <- xcell2object@signatures[startsWith(names(xcell2object@signatures), paste0(ctoi, "#"))]
    scores <- scoreMix(ctoi, mix_ranked, signatures_ctoi)
    return(scores)

  }, BPPARAM = param)
  names(res_raw) <- sigs_celltypes

  res <- t(sapply(res_raw, function(ctoi){
    rowMeans(ctoi)
  }))


  if (raw_scores) {
    res <- round(res, 4)
    return(res)
  }else{
    # linear transformation
    res <- t(sapply(rownames(res), function(ctoi){

      ctoi_res <- res[ctoi,]

      # Linear transformation
      a <- pull(filter(xcell2object@params, celltype == ctoi), a)
      b <- pull(filter(xcell2object@params, celltype == ctoi), b)
      m <- pull(filter(xcell2object@params, celltype == ctoi), m)
      ctoi_res <- (ctoi_res^(1/b)) / a
      ctoi_res <- ctoi_res*m

      # Shift values
      ctoi_res <- ctoi_res - min(ctoi_res)

      return(ctoi_res)
    }))
    if (ref_is_sc) {
      warningCondition("Reference type is scRNA-Seq - Spillover correction is disabled.")
      res <- round(res, 3)
      return(res)
    }
  }


  if (spillover) {

    # Spillover correction
    spill_mat <- xcell2object@spill_mat * spillover_alpha
    diag(spill_mat) <- 1

    rows <- intersect(rownames(res), rownames(spill_mat))

    scores_corrected <- apply(res[rows, ], 2, function(x) pracma::lsqlincon(spill_mat[rows, rows], x, lb = 0))
    scores_corrected[scores_corrected < 0] <- 0
    #scores_corrected <- round(scores_corrected, 4)
    rownames(scores_corrected) <- rows

    scores_corrected <- round(scores_corrected, 3)
    return(scores_corrected)
  }else{
    res <- round(res, 3)
    return(res)
  }

}

