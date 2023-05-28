#' xCell2Analysis function
#'
#' This function takes a matrix of bulk gene expression data and a `xCell2Signatures` object containing a set of signatures as input. It performs downstream analysis to identify enriched cell types in the bulk sample.
#'
#' @import singscore
#' @import dplyr
#' @import tibble
#' @import GSEABase
#' @import purrr
#' @param bulk A matrix containing gene expression data.
#' @param xcell2sigs A `xCell2Signatures` object.
#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
#' @examples
#' @export
xCell2Analysis <- function(bulk, xcell2sigs){

  # score ranked bulk gene expression matrix
  scoreBulk <- function(ctoi, bulk_ranked, xcell2sigs, genes_overlap = 0.8){

    signatures_ctoi <- xcell2sigs@filtered_signatures[startsWith(names(xcell2sigs@filtered_signatures), paste0(ctoi, "#"))]

    # Filter signatures by genes_overlap
    sig2remove <- c()
    for (i in 1:length(signatures_ctoi)) {
      sig <- signatures_ctoi[[i]]
      sig_genes <- GSEABase::geneIds(sig)
      shared_genes <- intersect(sig_genes, rownames(bulk_ranked))

      if (length(shared_genes) / length(sig_genes) < genes_overlap) {
        sig2remove <- c(sig2remove, names(signatures_ctoi)[i])
      }
    }

    if (length(sig2remove) == length(signatures_ctoi)) {
      warning(paste0("Not enough shared genes between bulk sample and ",  ctoi, " signatures - removing cell type..."))
      return(NA)
    }


    signatures_ctoi <- signatures_ctoi[!names(signatures_ctoi) %in% sig2remove]

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })

    median_scores <- apply(scores, 1, median)
    names(median_scores) <- colnames(bulk_ranked)

    return(median_scores)
  }

  # Transform scores
  transfomScores <- function(ctoi, scores, xcell2sigs){

    shift_value <- filter(xcell2sigs@transformation_parameters, celltype == ctoi)$shift_value
    shifted_scores <- scores - shift_value

    scaling_value <- filter(xcell2sigs@transformation_parameters, celltype == ctoi)$scaling_value
    transformation_model <- filter(xcell2sigs@transformation_parameters, celltype == ctoi)$model[[1]]
    transformed_scores <- round(predict(transformation_model, newdata = data.frame("shifted_score" = shifted_scores)), 4) * scaling_value

    return(transformed_scores)
  }

  # Rank bulk gene expression matrix
  bulk_ranked <- singscore::rankGenes(bulk)

  xCell2_out.tbl <- xcell2sigs@labels %>%
    as_tibble() %>%
    select(label = 2) %>%
    unique() %>%
    rowwise() %>%
    mutate(scores = list(scoreBulk(ctoi = label, bulk_ranked, xcell2sigs))) %>%
    filter(all(!is.na(scores))) %>%
    mutate(transfomed_scores = list(transfomScores(ctoi = label, scores, xcell2sigs))) %>%
    unnest_longer(transfomed_scores, indices_to = "sample") %>%
    select(-scores) %>%
    pivot_wider(names_from = sample, values_from = transfomed_scores)

  # Convert to matrix
  xCell2_out.mat <- as.matrix(xCell2_out.tbl[,-1])
  row.names(xCell2_out.mat) <- pull(xCell2_out.tbl[,1])

  return(xCell2_out.mat)
}


