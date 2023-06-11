# Signatures filtering analysis for scRNA-seq benchmarking

library(tidyverse)
library(xCell2)
data("ts_labels_with_ontology")

# Load dataset
datasets <- sub(".rds", "", list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/"))
ds <- "Thymus"

ds.rds <- readRDS(paste0("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/", ds, ".rds"))
ref <- ds.rds$singleCellExpr
labels <- tibble(ont = "NA",
                   label = ds.rds$singleCellLabels,
                   sample = colnames(ref),
                   dataset = ds.rds$singleCellSubjects) %>%
  mutate(ont = plyr::mapvalues(sample, ts_labels_with_ontology$sample, ts_labels_with_ontology$ont.fine, warn_missing = FALSE),
         label =  gsub("_", "-", label))
labels <- as.data.frame(labels)
bulk <- ds.rds$bulk
truth <- ds.rds$bulkRatio


# Generate signatures for dataset:
# (!) Run first functions of xCell2

# Get signatures correlations
getSigsCors <- function(sigs, bulk, truth, labels, method_info){

  bulk_ranked <- singscore::rankGenes(bulk)
  ds_celltypes <- unique(labels$label)


  # Score signatures
  ds_scores_list <- pbapply::pbsapply(ds_celltypes, function(ctoi){
    signatures_ctoi <- sigs[startsWith(names(sigs), paste0(ctoi, "#"))]
    scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    }))
    colnames(scores) <- colnames(bulk_ranked)
    rownames(scores) <- names(signatures_ctoi)
    scores
  })

  # Get correlations
  rownames(truth) <- gsub("_", "-", rownames(truth))
  ds_cors_list <- pbapply::pbsapply(names(ds_scores_list), function(ctoi){
    ds_ctoi_sigs_scores <- ds_scores_list[[ctoi]]

    sigs_cors <- apply(ds_ctoi_sigs_scores, 1, function(sig_scores){
      cor(sig_scores, truth[ctoi,], method = "spearman")
    })
  })

  # Make tidy
  enframe(ds_cors_list, name = "celltype") %>%
    unnest_longer(value, indices_to = "signature", values_to = "cor") %>%
    mutate(method = method_info) %>%
    return(.)
}
sigs_cors <- getSigsCors(sigs = signatures_collection, bulk, truth, labels, method_info = "none")


# Plot correlations
sigs_cors %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(position=position_jitter(width=0.2), alpha = 1/5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5)

sigs_cors %>%
  top_frac(n = 0.8, wt = cor) %>%
  group_by(celltype) %>%
  summarise(cor = mean(cor)) %>%
  summarise(mean_cor = mean(cor))




# Make mixture matrix (from pure_ct_mat) for first filtering creteria
mix_frac <- mixture_fractions[length(mixture_fractions)] # Use the highest fraction
makeMixutre <- function(pure_ct_mat, cor_mat, dep_list, mix_frac){

  # Make fractional CTOI and control matrices
  ctoi_mat <- pure_ct_mat * mix_frac
  celltypes <- colnames(pure_ct_mat)

  controls <- unname(sapply(celltypes, function(ctoi){
    dep_cts <- unname(unlist(dep_list[[ctoi]]))
    not_dep_cts <- celltypes[!celltypes %in% dep_cts]
    names(sort(cor_mat[ctoi, not_dep_cts])[1])
    }))

  controls_mat <- sapply(controls, function(ctrl){
    pure_ct_mat[,ctrl] * (1-mix_frac)
  })

  # Combine fractional matrices to a mixture
  mix_names <- paste0(colnames(ctoi_mat), "%%", colnames(controls_mat))
  mix_mat <- ctoi_mat + controls_mat
  colnames(mix_mat) <- mix_names

  # In case there is one control for all other cell types -> make a second controls matrix just for him
  if(length(unique(controls)) == 2){
    abundant_control <- names(sort(table(controls))[2])

    controls2 <- unname(sapply(celltypes, function(ctoi){
      dep_cts <- unname(unlist(dep_list[[ctoi]]))
      not_dep_cts <- celltypes[!celltypes %in% dep_cts]
      not_dep_cts <- not_dep_cts[not_dep_cts != abundant_control]
      names(sort(cor_mat[ctoi, not_dep_cts])[1])
    }))
    controls_mat2 <- sapply(controls2, function(ctrl){
      pure_ct_mat[,ctrl] * (1-mix_frac)
    })

    mix_names2 <- paste0(colnames(ctoi_mat), "%%", colnames(controls_mat2))
    mix_mat2 <- ctoi_mat + controls_mat2
    colnames(mix_mat2) <- mix_names2

    mixtures_list <- list(mix1 = mix_mat, mix2 = mix_mat2)

  }else{
    mixtures_list <- list(mix1 = mix_mat, mix2 = NULL)

  }

  return(mixtures_list)

}

mix_list <- makeMixutre(pure_ct_mat, cor_mat, dep_list, mix_frac)

scoreMixture <- function(mix, signatures_collection, dep_list){

  # Score
  mix_ranked <- singscore::rankGenes(mix)
  scores <- sapply(signatures_collection, simplify = TRUE, function(sig){
    singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })
  rownames(scores) <- colnames(mix)
  colnames(scores) <- names(signatures_collection)

  # Make tidy
  scores_tidy <- as_tibble(scores, rownames = "cts_sims") %>%
    pivot_longer(cols = -cts_sims, names_to = "signature", values_to = "score") %>%
    separate(cts_sims, into = c("mix_celltype", "mix_control"), sep = "%%") %>%
    separate(signature, into = "sig_celltype", sep = "#", extra = "drop", remove = FALSE)

  # Clean scores
  scores_tidy_clean <- scores_tidy %>%
    filter(sig_celltype != mix_control) %>% # Signature cell type cannot be the same as the control
    rowwise() %>%
    filter(!mix_celltype %in% unname(unlist(dep_list[[sig_celltype]])) & !mix_control %in% unname(unlist(dep_list[[sig_celltype]])))  # Simulation CTOI/control cannot be dependent on signature cell type

  return(scores_tidy_clean)

}




# First filtering - Top of delta score between CTOI and median score (of all other cell types)
scores_tidy_clean <- scoreMixture(mix_list$mix1, signatures_collection, dep_list)

median_sig_score <- scores_tidy_clean %>%
  ungroup() %>%
  filter(sig_celltype != mix_celltype) %>%
  group_by(signature) %>%
  summarise(median_score = median(score))

sigsPassed <- scores_tidy_clean %>%
  ungroup() %>%
  filter(sig_celltype == mix_celltype) %>%
  left_join(median_sig_score, by = "signature") %>%
  drop_na() %>%  # Some CTOI might not have scores with other cell types (because of controls)
  ungroup() %>%
  mutate(delta_score = score - median_score) %>%
  group_by(sig_celltype) %>%
  top_frac(n=0.2, wt=delta_score) %>%
  pull(signature)

# If cell types lost in this filtering step because of controls
if(!is.null(mix_list$mix2)){

  celltypes <- unique(labels$label)
  lost_ct <- celltypes[!celltypes %in% gsub("#.*", "", sigsPassed)]

  scores_tidy_clean_lost_ct <- scoreMixture(mix_list$mix2, signatures_collection[startsWith(names(signatures_collection), paste0(lost_ct, "#"))], dep_list)


  median_sig_score_lost_ct <- scores_tidy_clean_lost_ct %>%
    ungroup() %>%
    filter(sig_celltype != mix_celltype) %>%
    group_by(signature) %>%
    summarise(median_score = median(score))

  sigsPassed_lost_ct <- scores_tidy_clean_lost_ct %>%
    ungroup() %>%
    filter(sig_celltype == mix_celltype) %>%
    left_join(median_sig_score_lost_ct, by = "signature") %>%
    mutate(delta_score = score - median_score) %>%
    group_by(sig_celltype) %>%
    top_frac(n=0.2, wt=delta_score) %>%
    pull(signature)

  sigsPassed <- c(sigsPassed, sigsPassed_lost_ct)
}




sigs_cors %>%
  mutate(passed_filter = factor(ifelse(signature %in% sigsPassed , "yes", "no"), levels = c("yes", "no"))) %>%
  ggplot(., aes(x=celltype, y=cor, fill=passed_filter)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#66CD00", "#CD3333")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "green", size=0.5) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.2))

sigs_cors %>%
  filter(signature %in% sigsPassed) %>%
  group_by(celltype) %>%
  summarise(cor = mean(cor)) %>%
  summarise(mean_cor = mean(cor))


# Second filtering -
signatures_filtered <- signatures_collection[names(signatures_collection) %in% sigsPassed]

# Make simulations
makeSimulations <- function(ref, labels, mixture_fractions, dep_list, n_ct_sim = 10, add_noise = FALSE, seed = 123){

  set.seed(seed)

  makeFractionMatrixCTOI <- function(ref, mixture_fractions, ctoi_samples2use){


    ctoi_median_expression <- Rfast::rowMedians(ref[,ctoi_samples2use])
    ctoi_frac_mat <- matrix(rep(ctoi_median_expression, length(mixture_fractions)), byrow = FALSE, ncol = length(mixture_fractions)) %*% diag(mixture_fractions)
    rownames(ctoi_frac_mat) <- rownames(ref)

    # # Make CTOI fraction matrix
    # ctoi_frac_mat <- matrix(NA, nrow = nrow(ref), ncol = length(mixture_fractions), dimnames = list(rownames(ref), mixture_fractions))
    # for (i in 1:length(mixture_fractions)) {
    #   frac <- mixture_fractions[i]
    #   frac_fracs <- diff(c(0, sort(runif(length(ctoi_samples2use)-1, min = 0, max = frac)), frac)) # Generate random fraction for each sample to sum to frac in mixture_fractions
    #   ctoi_frac_mat[,i] <- Rfast::rowsums(ref[,ctoi_samples2use] %*% diag(frac_fracs)) # Sum all fractions
    # }

    return(ctoi_frac_mat)
  }
  makeFractionMatrixControls <- function(ref, labels, mixture_fractions, control_cts){


    # Pick control samples
    control_samples <- labels %>%
      filter(label == control_cts) %>%
      slice_sample(n=length(mixture_fractions)) %>%
      pull(sample)

    if (length(control_samples) == 1) {
      controls_median_expression <- ref[,control_samples]
    }else{
      controls_median_expression <- Rfast::rowMedians(ref[,control_samples])
    }

    controls_mat <- matrix(rep(controls_median_expression, length(mixture_fractions)), byrow = FALSE, ncol = length(mixture_fractions)) %*% diag(1-mixture_fractions)
    rownames(controls_mat) <- rownames(ref)


    # controls_mat <- sapply(1-mixture_fractions, function(frac){
    #   random_fracs <- diff(c(0, sort(runif(ncol(controls_expression)-1, min = 0, max = frac)), frac)) # Generate random numbers from the controls from a uniform distribution that sum to frac
    #   Rfast::rowsums(controls_expression %*% diag(random_fracs))
    # })

    return(controls_mat)

  }


  celltypes <- unique(labels[,2])
  sim_list <- pbapply::pblapply(celltypes, function(ctoi){


    # Sort CTOI samples to be homogeneous by datasets
    ctoi_samples_pool <- c()
    while(!all(labels[labels$label == ctoi,]$sample %in% ctoi_samples_pool)) {
      ctoi_samples_pool <- c(ctoi_samples_pool,
                             labels %>%
                               filter(label == ctoi & !sample %in% ctoi_samples_pool) %>%
                               slice_head(n = 1, by = dataset) %>%
                               pull(sample))
    }

    if (length(ctoi_samples_pool) < length(mixture_fractions)) {
      ctoi_samples_pool <- rep(ctoi_samples_pool, length(mixture_fractions))
    }

    # Get dependent cell types
    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))

    # Get number of available datasets of this cell type (for adding noise)...
    n_ctoi_ds <- length(unique(labels[labels$label == ctoi,]$dataset))


    ctoi_sim_list <- lapply(1:n_ct_sim, function(i){

      ctoi_samples2use <- ctoi_samples_pool[1:length(mixture_fractions)] # Choose the first samples (homogeneous by datasets)
      ctoi_frac_mat <- makeFractionMatrixCTOI(ref, mixture_fractions, ctoi_samples2use)
      ctoi_samples_pool <- c(ctoi_samples_pool[!ctoi_samples_pool %in% ctoi_samples2use], ctoi_samples2use) # Move ctoi_samples2use to be last

      # Make Controls fraction matrix
      control_ct <- names(sort(cor_mat[ctoi, !colnames(cor_mat) %in% dep_cts])[1])
      controls_frac_mat <- makeFractionMatrixControls(ref, labels, mixture_fractions, control_ct)

      # Combine CTOI and controls fractions matrix
      simulation <- ctoi_frac_mat + controls_frac_mat
      colnames(simulation) <- paste0(control_ct, "%%", mixture_fractions)

      # Add noise
      if (add_noise) {
        noise_sd <- 1/n_ctoi_ds
        noise <- matrix(rnorm(nrow(simulation) * ncol(simulation), mean = 0, sd = noise_sd),
                        nrow = nrow(simulation), ncol = ncol(simulation))
        simulation <- simulation + noise
        simulation <- pmax(simulation, 0)
      }

      simulation
    })
    names(ctoi_sim_list) <- paste0("sim-", 1:n_ct_sim)

    ctoi_sim_list

  })
  names(sim_list) <- celltypes

  return(sim_list)

}
sim_list <- makeSimulations(ref, labels, mixture_fractions, dep_list, n_ct_sim = 10, add_noise = TRUE, seed = 123)

# Score CTOI simulations
scoreCTOISimulations <- function(signatures, sim_list){

  celltypes <- names(sim_list)
  out <- pbapply::pblapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

    scores_ctoi <- lapply(sim_list[[ctoi]], function(sim){
      sim_ranked <- singscore::rankGenes(sim)
      scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      }))

      rownames(scores) <- names(signatures_ctoi)
      colnames(scores) <- colnames(sim)

      scores
    })
    scores_ctoi

  })
  names(out) <- celltypes
  return(out)
}
scores_list <- scoreCTOISimulations(signatures = signatures_filtered, sim_list)

# Get simulation correlations results
mat2tidy <- function(mat){
  names <- rownames(mat)
  as_tibble(mat) %>%
    mutate(signature=names) %>%
    relocate(signature) %>%
    pivot_longer(-signature, names_to = "fraction", values_to = "score") %>%
    #mutate(fraction = as.numeric(fraction)) %>%
    return(.)
}

scores_all_sims_tidy <- enframe(scores_list, name = "celltype") %>%
  unnest_longer(value, indices_to = "sim_id", values_to = "scores") %>%
  rowwise() %>%
  mutate(scores = list(mat2tidy(scores))) %>%
  unnest(scores) %>%
  separate(fraction, into = c("control", "fraction"), sep = "%%") %>%
  mutate(fraction = as.numeric(fraction))

# Filter signatures
sigsPassed2 <- scores_all_sims_tidy %>%
  group_by(celltype, signature, sim_id) %>%
  summarise(cor = cor(fraction, score, method = "spearman")) %>%
  group_by(celltype, signature) %>%
  summarise(mean_cor = mean(cor)) %>%
  top_frac(n=0.2, wt=mean_cor) %>%
  pull(signature)


sigs_cors %>%
  filter(signature %in% sigsPassed2) %>%
  group_by(celltype) %>%
  summarise(cor = mean(cor)) %>%
  summarise(mean_cor = mean(cor))


sigs_cors %>%
  mutate(passed_filter = factor(ifelse(signature %in% sigsPassed2 , "yes", "no"), levels = c("yes", "no"))) %>%
  ggplot(., aes(x=celltype, y=cor, fill=passed_filter)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#66CD00", "#CD3333")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "green", size=0.5) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.2))






