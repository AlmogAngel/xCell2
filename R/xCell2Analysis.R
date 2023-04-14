########################################################################################
# Analyze new mixtures
########################################################################################

library(tidyverse)

xCell2Analysis <- function(mix, xcell2sigs){

  scoreMixtures <- function(ctoi, mixture_ranked){
    signatures_ctoi <- xcell2sigs@filtered_signatures[startsWith(names(xcell2sigs@filtered_signatures), paste0(ctoi, "#"))]

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mixture_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })

    # In case some signatures contain genes that are all not in the mixtures
    if (is.list(scores)) {
      signatures_ctoi <- signatures_ctoi[-which(lengths(scores) == 0)]
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(mixture_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      })
    }

    rownames(scores) <- colnames(mixture_ranked)
    return(scores)
  }


  # Rank mixture genes
  mix_ranked <- singscore::rankGenes(mix)

  # By mean scores
  xCell2_out <- xcell2sigs@labels %>%
    as_tibble() %>%
    select(label = 2) %>%
    unique() %>%
    rowwise() %>%
    mutate(scores = list(rowMeans(scoreMixtures(ctoi = label, mix_ranked)))) %>%
    mutate(samples = list(names(scores))) %>%
    unnest(cols = c(samples, scores)) %>%
    pivot_wider(names_from = samples, values_from = scores) %>%
    as.data.frame()


  # # By models
  # xCell2_out <- xcell2sigs@models %>%
  #   rowwise() %>%
  #   mutate(scores = list(scoreMixtures(ctoi = label, mixture = mix_ranked, signatures_ctoi = signatures))) %>%
  #   # For lasso:
  #   #mutate(predictions = list(as.numeric(predict(lasso, newx = scores)))) %>%
  #   # For GGRF:
  #   #mutate(predictions = list(as.numeric(predict(GRRF, newdata = scores)))) %>%
  #   dplyr::select(label, predictions) %>%
  #   unnest(predictions) %>%
  #   mutate(samples = rep(colnames(mix), length(celltypes))) %>%
  #   pivot_wider(names_from = samples, values_from = predictions) %>%
  #   as.data.frame()

  rownames(xCell2_out) <- xCell2_out[,1]
  xCell2_out <- xCell2_out[,-1]

  return(xCell2_out)
}
