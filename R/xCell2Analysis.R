#' xCell2Analysis function
#'
#' This function takes a matrix of mix gene expression data and a `xCell2Signatures` object containing a set of signatures as input. It performs downstream analysis to identify enriched cell types in the mix sample.
#'
#' @importFrom singscore rankGenes simpleScore
#' @import dplyr
#' @import tibble
#' @import purrr
#' @importFrom  pracma lsqlincon
#' @param mix a matrix containing gene expression data
#' @param xcell2sigs S4 object of `xCell2Signatures`
#' @param tranform Boolean - should scores transformed to fraction?
#' @param spillover Boolean - should we use spillover corretion on the transformed scores?
#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
#' @examples
#' @export
xCell2Analysis <- function(mix, xcell2sigs, min_intersect = 0.9, tranform, spillover, spillover_alpha = 0.5){

  # score ranked mix gene expression matrix
  scoreMix <- function(ctoi, mix_ranked, xcell2sigs){

    signatures_ctoi <- xcell2sigs@signatures[startsWith(names(xcell2sigs@signatures), paste0(ctoi, "#"))]

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
    })

    return(scores)
  }

  predictFracs <- function(ctoi, scores, xcell2sigs){

    model <- filter(xcell2sigs@transformation_models, celltype == ctoi)$model[[1]]
    predictions <- round(predict(model, newdata = scores, type = "response"), 3)

    return(predictions)

  }

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

   all_predictions <- tibble(label = sigs_celltypes) %>%
    rowwise() %>%
    mutate(scores = list(scoreMix(ctoi = label, mix_ranked, xcell2sigs))) %>%
    mutate(predictions = ifelse(tranform, list(predictFracs(ctoi = label, scores, xcell2sigs)), list(scores))) %>%
    unnest_longer(predictions, indices_to = "sample") %>%
    select(-scores) %>%
    pivot_wider(names_from = sample, values_from = predictions)

  # Convert to matrix
   all_predictions_mat <- as.matrix(all_predictions[,-1])
  rownames(all_predictions_mat) <- pull(all_predictions[,1])

  if (tranform & spillover) {
    # Spillover correction
    spill_mat <- xcell2sigs@spill_mat * spillover_alpha
    diag(spill_mat) <- 1

    rows <- rownames(scores_transformed_mat)[rownames(scores_transformed_mat) %in% rownames(spill_mat)]

    scores_corrected <- apply(scores_transformed_mat[rows, ], 2, function(x) pracma::lsqlincon(spill_mat[rows, rows], x, lb = 0))
    scores_corrected[scores_corrected < 0] <- 0
    rownames(scores_corrected) <- rows
    return(scores_corrected)
  }else{
    return(scores_transformed_mat)
  }



}


