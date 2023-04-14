library(tidyverse)

getMixtures <- function(ctoi, ref_ctoi, pure_ct_mat, dep.list, fracs){

  fractions <- sort(rep(fracs, ncol(ref_ctoi)))

  # Make a matrix with all fractions for all cell type samples
  ref_ctoi.rep <- matrix(rep(ref_ctoi, length(fracs)), ncol = length(fracs)*ncol(ref_ctoi),
                         byrow = FALSE, dimnames = list(rownames(ref_ctoi), paste0(rep(colnames(ref_ctoi), length(fracs)), "#", fractions))) # Replicated matrix by length of fractions
  ref_ctoi.rep.frac <- t(t(ref_ctoi.rep) * fractions) # Multiply matrix by fractions

  # add background cell type expression
  if (length(dep.list[[ctoi]]) == 0) {
    not_dep_cells <- rep(TRUE, ncol(pure_ct_mat))
  }else{
    not_dep_cells <- !colnames(pure_ct_mat) %in% dep.list[[ctoi]]
  }

  if (sum(not_dep_cells) > 2) {
    median_back_exp <- Rfast::rowMedians(pure_ct_mat[,not_dep_cells])
  }else{
    median_back_exp <- pure_ct_mat[,not_dep_cells]
  }


  median_back_exp.mat <- matrix(rep(median_back_exp, ncol(ref_ctoi.rep.frac)), ncol = ncol(ref_ctoi.rep.frac), byrow = FALSE)
  median_back_exp.mat.frac <- t(t(median_back_exp.mat) * (1-fractions)) # Multiply matrix by 1 - fractions

  ctoi_mix <- ref_ctoi.rep.frac + median_back_exp.mat.frac

  return(list(mixture = ctoi_mix, fractions = fractions))
}


scoreMixtures <- function(ctoi, mixture, signatures_ctoi){

  mixture_ranked <-  singscore::rankGenes(mixture)
  scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
    singscore::simpleScore(mixture_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })


  return(scores)
}


trainModels <- function(ref, labels, dep_list, pure_ct_mat_test, signatures_collection_filtered, mixture_fractions){

  # Build mixtures
  mixtures <- labels %>%
    as_tibble() %>%
    dplyr::select(label) %>%
    unique() %>%
    rowwise() %>%
    mutate(ref_label = list(ref[,labels$is_test & labels$label == label])) %>%
    mutate(mixtures = list(getMixtures(ctoi = label, ref_ctoi = ref_label, pure_ct_mat = pure_ct_mat_test, dep.list = dep_list, fracs = mixture_fractions)))


  # Score mixtures
  mixtures.scored <- mixtures %>%
    rowwise() %>%
    mutate(signatures = list(signatures_collection_filtered[startsWith(names(signatures_collection_filtered), paste0(label, "#"))])) %>%
    mutate(scores = list(scoreMixtures(ctoi = label, mixture = mixtures$mixture, signatures_ctoi = signatures)))


  # # Train models (lasso)
  # mixtures.models <- mixtures.scored %>%
  #   rowwise() %>%
  #   mutate(lasso = list(glmnet::cv.glmnet(scores, mixtures$fractions, alpha = 1, family = "gaussian"))) %>%
  #   dplyr::select(label, signatures, lasso)

  # Train models (GRRF) https://sites.google.com/site/houtaodeng/rrf
  gamma <- 0.8 # A larger gamma often leads to fewer features. But, the quality of the features selected is quite stable for GRRF
  mixtures.models <- mixtures.scored %>%
    rowwise() %>%
    mutate(RF = list(RRF::RRF(scores, mixtures$fractions, flagReg = 0, importance = TRUE))) %>%
    mutate(RF_imp = list(RF$importance[,"%IncMSE"] / max(RF$importance[,"%IncMSE"]))) %>%
    mutate(GRRF = list(RRF::RRF(scores, mixtures$fractions, flagReg = 1, coefReg = (1-gamma) + gamma*RF_imp))) %>%
    dplyr::select(label, signatures, GRRF)

  return(mixtures.models)
}




#mixtures.models[1,]$signatures
#coef(mixtures.models[1,]$lasso[[1]])












# # Add background cell types to the pure CTOI matrix
# getFracs <- function(n_ct_mix, ctoi_frac){
#   x <- runif(n_ct_mix, 0, 1)
#   fracs <- x * (1-ctoi_frac) / sum(x)
#   return(round(fracs, 4))
# }


# ref_ctoi <- mixtures[1,]$ref_ctoi[[1]]
#
# getMixtures <- function(ctoi, ref_ctoi, mixture_fractions, dep_list){
#   pure_ct_mat.nodep <- pure_ct_mat_train[,!colnames(pure_ct_mat_train) %in% dep_list[[ctoi]]] # pure_ct_mat_train with no dependencies
#
#   for (i in 1:ncol(ref_ctoi)) {
#     # For each sample
#     sample <- ref_ctoi[,i]
#     ref_ctoi.fracs <- matrix(rep(sample, length(mixture_fractions)), ncol = length(mixture_fractions), byrow = FALSE) # Replicated matrix by length of fractions
#     ref_ctoi.fracs <- t(t(ref_ctoi.fracs) * mixture_fractions) # Multiply matrix by fractions
#
#     # add background cell type expression
#     ref_ctoi.fracs.reps <- ref_ctoi.fracs
#     for (rep in 1:n_reps) {
#       random_type_samples <- sample(1:ncol(pure_ct_mat.nodep), n_ct_mix) # Get n_ct_mix random background cell types for the mixtures
#       ref_ctoi.fracs.rep <- ref_ctoi.fracs
#
#       for (i in 1:ncol(ref_ctoi.fracs)) {
#         remaining_fracs <-  getFracs(n_ct_mix, mixture_fractions[i]) # Get n_ct_mix random fractions to complete ctoi fraction
#         background_median <- Rfast::rowMedians(t(t(pure_ct_mat.nodep[,random_type_samples]) * remaining_fracs)) # Get the median expression of the fraction of the backgroud cell-types
#         ref_ctoi.fracs.rep[,i] <- ref_ctoi.fracs.rep[,i] + background_median # Add to the ctoi fraction
#       }
#
#       ref_ctoi.fracs.reps <- cbind(ref_ctoi.fracs.reps, ref_ctoi.fracs.rep)
#     }
#     ref_ctoi.fracs.reps <- ref_ctoi.fracs.reps[,-c(1:ncol(ref_ctoi.fracs))]
#
#
#
#   }
#
# }











