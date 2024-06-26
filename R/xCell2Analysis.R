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
#' @param min_intersect description
#' @param raw_scores description
#' @param estimate_fracs Boolean
#' @param spillover Boolean - should we use spillover correction on the transformed scores?
#' @param spillover_alpha description
#' @param num_threads description

#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
#' @export
xCell2Analysis <- function(mix, xcell2object, min_intersect = 0.9, raw_scores = FALSE, estimate_fracs = FALSE, spillover = TRUE, spillover_alpha = 0.2, num_threads = 1){

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
  if (genes_intersect_frac < min_intersect) {
    stop("Intersect between reference and mixture's genes is: ", genes_intersect_frac, " and min_intersect = ", min_intersect, ".",
         "\n", "Please use xCell2CleanGenes before using xCell2 signatures with that mixture.")
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

  if (estimate_fracs) {

    # a <- pull(filter(xcell2object@params, celltype == ctoi), a)
    # b <- pull(filter(xcell2object@params, celltype == ctoi), b)
    # intercepts <- pull(filter(xcell2object@params, celltype == ctoi), intercepts)[[1]]
    # betas <- pull(filter(xcell2object@params, celltype == ctoi), betas)[[1]]
    #
    # # Linear transformation
    # scores <- (scores^(1/b)) / a
    #
    # # Scale
    # scores <- scale(scores)
    #
    # # Predict
    # p <- (scores %*% betas) + intercepts
    # p <- apply(p, 1, mean)
    #
    #
    # return(p)
  }else{
    res <- t(sapply(res_raw, function(ctoi){
      rowMeans(ctoi)
    }))

  }

  if (raw_scores) {
    res <- round(res, 3)
    return(res)
  }else{
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

      ctoi_res <- round(ctoi_res, 5)
      return(ctoi_res)
    }))
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

    return(scores_corrected)
  }else{
    return(res)
  }

}

