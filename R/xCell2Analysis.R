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
#' @param xcell2sigs S4 object of `xCell2Signatures`
#' @param predict Boolean - should we use model for to predict final scores?
#' @param spillover Boolean - should we use spillover corretion on the transformed scores?
#' @param ncores

#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
#' @export
xCell2Analysis <- function(mix, xcell2sigs, min_intersect = 0.9, predict, spillover, spillover_alpha = 0.2, ncores = 1){

  scoreMix <- function(ctoi, mix_ranked, signatures_ctoi){

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
    })

    # Shift values
    scores <- apply(scores, 2, function(s){
      s - min(s)
    })

    rownames(scores) <- colnames(mix_ranked)

    return(scores)
  }
  predictFracs <- function(ctoi, scores, model_ctoi){

    # predictions <- round(randomForestSRC::predict.rfsrc(model, newdata = as.data.frame(scores))$predicted, 4)
    # predictions <- round(predict(model, scores, s = model$lambda, type = "response")[,1], 4)
    predictions <- round(predict(model_ctoi, scores, type = "response"), 8)


    names(predictions) <- rownames(scores)

    return(predictions)

  }

  param <- BiocParallel::MulticoreParam(workers = ncores)


  # Check reference/mixture genes intersection
  genes_intersect_frac <- length(intersect(rownames(mix), xcell2sigs@genes_used)) / length(xcell2sigs@genes_used)
  if (genes_intersect_frac < min_intersect) {
    stop("Intersect between reference and mixture's genes is: ", genes_intersect_frac, " and min_intersect = ", min_intersect, ".",
         "\n", "Please use xCell2CleanGenes before using xCell2 signatures with that mixture.")
  }


  # Rank mix gene expression matrix
  mix_ranked <- singscore::rankGenes(mix)

  # Score and predict
  sigs_celltypes <- unique(unlist(lapply(names(xcell2sigs@signatures), function(x){strsplit(x, "#")[[1]][1]})))

  res <- BiocParallel::bplapply(sigs_celltypes, function(ctoi){

    signatures_ctoi <- xcell2sigs@signatures[startsWith(names(xcell2sigs@signatures), paste0(ctoi, "#"))]
    model_ctoi <- filter(xcell2sigs@models, celltype == ctoi)$model[[1]]

    scores <- scoreMix(ctoi, mix_ranked, signatures_ctoi)
    if (predict) {
      return(predictFracs(ctoi, scores, model_ctoi))
    }else{
      return(rowMeans(scores))
    }

  }, BPPARAM = param)
  names(res) <- sigs_celltypes

  res <- t(sapply(res, cbind))
  colnames(res) <-  colnames(mix_ranked)


  if (!spillover) {
    return(res)
  }


  # Spillover correction
  spill_mat <- xcell2sigs@spill_mat * spillover_alpha
  diag(spill_mat) <- 1

  rows <- intersect(rownames(res), rownames(spill_mat))

  scores_corrected <- apply(res[rows, ], 2, function(x) pracma::lsqlincon(spill_mat[rows, rows], x, lb = 0))
  scores_corrected[scores_corrected < 0] <- 0
  scores_corrected <- round(scores_corrected, 4)
  rownames(scores_corrected) <- rows
  return(scores_corrected)


}

