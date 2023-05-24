#' xCell2Analysis function
#'
#' This function takes a matrix of bulk gene expression data and a `xCell2Signatures` object containing a set of signatures as input. It performs downstream analysis to identify enriched cell types in the bulk sample.
#'
#' @import singscore
#' @import dplyr
#' @import tibble
#' @import GSEABase
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

    scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    }))


    colnames(scores) <- colnames(bulk_ranked)
    rownames(scores) <- names(signatures_ctoi)

    return(scores)
  }

  # Transform scores
  transfomScores <- function(ctoi, scores, xcell2sigs){

    xcell2sigs@transformation_parameters[xcell2sigs@transformation_parameters$celltype == ctoi,] %>%
      rowwise() %>%
      mutate(scores = if(signature %in% rownames(scores)) list(scores[signature,]) else NA) %>%
      print()
      # drop_na() %>%
      # mutate(shifted_score = map2(shift_value, scores, ~ .y - .x)) %>%  # Shift scores
      # mutate(shifted_score = map(shifted_score, ~ pmax(., 0))) %>%
      # rowwise() %>%
      # mutate(trasfomed_score = list(round(predict(betareg, newdata = data.frame(shifted_score), type = "response")*scaling_value, 3))) %>%
      # unnest_longer(trasfomed_score, indices_to = "sample") %>%
      # group_by(celltype, sample) %>%
      # summarise(fraction = mean(trasfomed_score)) %>%
      # ungroup() %>%
      # select(-celltype) %>%
      return(.)
  }

  # Rank bulk gene expression matrix
  bulk_ranked <- singscore::rankGenes(bulk)

  xCell2_out.tbl <- xcell2sigs@labels %>%
    as_tibble() %>%
    select(label = 2) %>%
    unique() %>%
    rowwise() %>%
    mutate(scores = list(scoreBulk(ctoi = label, bulk_ranked, xcell2sigs))) %>%
    filter(all(!is.na(scores))) %>%  # Remove cell types with low gene overlap with bulk sample
    mutate(ct_fractions = list(transfomScores(ctoi = label, scores, xcell2sigs))) %>%
    unnest(ct_fractions) %>%
    select(-scores) %>%
    pivot_wider(names_from = sample, values_from = fraction)

  # Convert to matrix
  xCell2_out.mat <- as.matrix(xCell2_out.tbl[,-1])
  row.names(xCell2_out.mat) <- pull(xCell2_out.tbl[,1])

  return(xCell2_out.mat)
}
