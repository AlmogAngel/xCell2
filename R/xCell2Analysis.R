#' xCell2Analysis function
#'
#' This function takes a matrix of bulk gene expression data and a `xCell2Signatures` object containing a set of signatures as input. It performs downstream analysis to identify enriched cell types in the bulk sample.
#'
#' @importFrom singscore rankGenes simpleScore
#' @import dplyr
#' @import tibble
#' @import purrr
#' @importFrom  pracma lsqlincon
#' @param bulk A matrix containing gene expression data.
#' @param xcell2sigs A `xCell2Signatures` object.
#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
#' @examples
#' @export
xCell2Analysis <- function(bulk, xcell2sigs, min_genes_overlap = 0.5, spillover_alpha = 0.5){

  # score ranked bulk gene expression matrix
  scoreBulk <- function(ctoi, bulk_ranked, xcell2sigs, genes_overlap = min_genes_overlap){

    signatures_ctoi <- xcell2sigs@signatures[startsWith(names(xcell2sigs@signatures), paste0(ctoi, "#"))]

    # Filter signatures by genes_overlap
    sig2remove <- c()
    for (i in 1:length(signatures_ctoi)) {
      sig_genes <- signatures_ctoi[[i]]

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
      suppressWarnings(singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
    })

    mean_scores <- apply(scores, 1, mean)
    names(mean_scores) <- colnames(bulk_ranked)

    return(mean_scores)
  }

  # Transform scores
  transfomScores <- function(ctoi, scores, xcell2sigs){

    shift_value <- filter(xcell2sigs@transformation_models, celltype == ctoi)$shift_value
    shifted_scores <- scores - shift_value

    transformation_model <- filter(xcell2sigs@transformation_models, celltype == ctoi)$model[[1]]
    transformed_scores <- round(predict(transformation_model, newdata = data.frame("shifted_score" = shifted_scores), type = "response"), 2)

    return(transformed_scores)
  }

  # Rank bulk gene expression matrix
  bulk_ranked <- singscore::rankGenes(bulk)

  # Score and transform
  sigs_celltypes <- unique(unlist(lapply(names(xcell2sigs@signatures), function(x){strsplit(x, "#")[[1]][1]})))

  scores_transformed <- tibble(label = sigs_celltypes) %>%
    rowwise() %>%
    mutate(scores = list(scoreBulk(ctoi = label, bulk_ranked, xcell2sigs))) %>%
    filter(all(!is.na(scores))) %>%
    mutate(transfomed_scores = list(transfomScores(ctoi = label, scores, xcell2sigs))) %>%
    unnest_longer(transfomed_scores, indices_to = "sample") %>%
    select(-scores) %>%
    pivot_wider(names_from = sample, values_from = transfomed_scores)

  # Convert to matrix
  scores_transformed_mat <- as.matrix(scores_transformed[,-1])
  row.names(scores_transformed_mat) <- pull(scores_transformed[,1])

  # Spillover correction
  spill_mat <- xcell2sigs@spill_mat * spillover_alpha
  diag(spill_mat) <- 1

  rows <- rownames(scores_transformed_mat)[rownames(scores_transformed_mat) %in% rownames(spill_mat)]

  scores_corrected <- apply(scores_transformed_mat[rows, ], 2, function(x) pracma::lsqlincon(spill_mat[rows, rows], x, lb = 0))
  scores_corrected[scores_corrected < 0] <- 0
  rownames(scores_corrected) <- rows

  return(scores_corrected)
}


