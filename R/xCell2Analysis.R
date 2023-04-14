#' xCell2Analysis function
#'
#' This function takes a matrix of bulk gene expression data and a `xCell2 Signatures` object containing a set of signatures as input. It performs downstream analysis to identify enriched cell types in the bulk sample.
#'
#' @import singscore
#' @import dplyr
#' @import tibble
#' @param bulk A matrix containing gene expression data.
#' @param xcell2sigs A `xCell2 Signatures` object.
#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
#' @examples
xCell2Analysis <- function(bulk, xcell2sigs){

  # score ranked bulk gene expression matrix
  bulk_scored <- function(ctoi, bulk_ranked){
    signatures_ctoi <- xcell2sigs@filtered_signatures[startsWith(names(xcell2sigs@filtered_signatures), paste0(ctoi, "#"))]

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })

    # In case some signatures contain genes that are all not in the bulk matrix
    if (is.list(scores)) {
      signatures_ctoi <- signatures_ctoi[-which(lengths(scores) == 0)]
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      })
    }

    rownames(scores) <- colnames(bulk_ranked)
    return(scores)
  }

  # Rank bulk gene expression matrix
  bulk_ranked <- singscore::rankGenes(bulk)

  # By mean scores
  xCell2_out <- xcell2sigs@labels %>%
    as_tibble() %>%
    select(label = 2) %>%
    unique() %>%
    rowwise() %>%
    mutate(scores = list(rowMeans(bulk_scored(ctoi = label, bulk_ranked)))) %>%
    mutate(samples = list(names(scores))) %>%
    unnest(cols = c(samples, scores)) %>%
    pivot_wider(names_from = samples, values_from = scores) %>%
    as.data.frame()

  rownames(xCell2_out) <- xCell2_out[,1]
  xCell2_out <- xCell2_out[,-1]

  return(xCell2_out)
}
