library(tidyverse)
# Function for generating simulations

labels <- readRDS("/bigdata/almogangel/xCell2/dev_data/sref_blood_labels_bulk.rds")
ref <- readRDS("/bigdata/almogangel/xCell2/dev_data/sref_blood_data_bulk.rds")
dep_list <- xCell2::xCell2GetLineage(labels = labels[,1:2], out_file = NULL)
mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02))

makeSimulations <- function(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = FALSE){


  makeControlsMatrix <- function(ref, labels, ctoi, dep_cts, mixture_fractions, ds_used){

    # In case we need to take more samples from a dataset that already have been used
    if(all(labels[labels$label == ctoi,]$dataset %in% ds_used)){
      ds2use <- unique(labels$dataset)
    }else{
      ds2use <- labels[!labels$dataset %in% ds_used,]$dataset
    }

    # Pick controls
    control_samples <- labels %>%
      filter(!label %in% dep_cts) %>%
      filter(dataset %in% ds2use) %>%
      group_by(label) %>%
      slice_sample(n=1) %>%
      pull(sample)
    controls_expression <- ref[,control_samples]

    controls_mat <- sapply(1-mixture_fractions, function(frac){
      random_fracs <- diff(c(0, sort(runif(ncol(controls_expression)-1, min = 0, max = frac)), frac)) # Generate random numbers from the controls from a uniform distribution that sum to frac
      Rfast::rowsums(controls_expression %*% diag(random_fracs))
    })

    return(controls_mat)

  }

  # This function return a fraction matrix of pure cell type
  # Ideally from different samples in different datasets
  makeFractionMatrix <- function(ref, labels, ctoi, mixture_fractions, ds_used){


    if(all(labels[labels$label == ctoi,]$dataset %in% ds_used)){ # If all datasets have been used
      ds2use <- unique(labels$dataset)
    }else{
      ds2use <- labels[!labels$dataset %in% ds_used,]$dataset # Only datasets that haven't been used
    }

    # Pick samples
    samples_pool <- labels %>%
      filter(label == ctoi) %>%
      filter(dataset %in% ds2use)

    ctoi_samples <- c()
    if (nrow(samples_pool) < length(mixture_fractions)) {
      # If there are not enough samples from different datasets
      ctoi_samples <- rep(samples_pool$sample, length(mixture_fractions))[1:length(mixture_fractions)]
    }else{

      # This function pulls samples that havn't been used from different datasets
      getSamples <- function(samples_pool, samples_used){
        samples_pool %>%
          filter(!sample %in% samples_used) %>%
          slice_sample(n=1, by = dataset) %>%
          pull(sample) %>%
          return(.)
      }

      while (length(ctoi_samples) < length(mixture_fractions)) {
        ctoi_samples <- c(ctoi_samples, getSamples(samples_pool, samples_used = ctoi_samples))
      }
      ctoi_samples <- ctoi_samples[1:length(mixture_fractions)]
    }

    ctoi_expression_mat <- ref[,ctoi_samples]

    # Make fraction matrix
    frac_mat <- ctoi_expression_mat %*% diag(mixture_fractions)
    colnames(frac_mat) <- mixture_fractions
    rownames(frac_mat) <- rownames(ref)

    ds_used <- labels[labels$sample == ctoi_samples,]$dataset

    out <- list(frac_mat = frac_mat, ds_used = ds_used)

    return(out)
  }

  celltypes <- unique(labels[,2])

  # Generate a list of n_ct_sim fraction matrices for each cell type
  sim_list <- pbapply::pblapply(celltypes, function(ctoi){
    print(ctoi)

    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
    ctoi_sim_list <- list()
    n_ctoi_ds <- length(unique(labels[labels$label == ctoi,]$dataset))
    ds_used <- c()
    for (i in 1:n_ct_sim) {

      # Make CTOI fractions matrix
      makeFractionMatrix.out <- makeFractionMatrix(ref, labels, ctoi, mixture_fractions, ds_used)
      ctoi_frac_mat <- makeFractionMatrix.out$frac_mat

      # Make Controls 1-fraction matrix
      controls_frac_mat <- makeControlsMatrix(ref, labels, ctoi, dep_cts, mixture_fractions, ds_used)

      # Combine CTOI and controls fractions matrix
      simulation <- ctoi_frac_mat + controls_frac_mat

      # Add noise
      if (add_noise & n_ctoi_ds < 10) {

        if (n_ctoi_ds == 1) {
          noise_sd <- 1.5
        }else if(n_ctoi_ds > 1 & n_ctoi_ds < 5){
          noise_sd <- 1
        }else if(n_ctoi_ds >= 5){
          noise_sd <- 0.5
        }

        noise <- matrix(rnorm(nrow(simulation) * ncol(simulation), mean = 0, sd = noise_sd),
                        nrow = nrow(simulation), ncol = ncol(simulation))
        simulation <- combined_mat + noise
      }

      ctoi_sim_list[[paste0("sim-", i)]] <- simulation
      ds_used <- c(ds_used, makeFractionMatrix.out$ds_used)
    }

    ctoi_sim_list
  })
  names(sim_list) <- celltypes

  return(sim_list)

}

# Make simulation with and without noise

sim_list <- makeSimulations(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = FALSE)
# saveRDS(sim_list, "/bigdata/almogangel/xCell2/dev_data/sim_list080523.rds")

sim_list_noise_sd0.5 <- makeSimulations(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = TRUE, noise_sd = 0.5)
# saveRDS(sim_list_noise_sd0.5, "/bigdata/almogangel/xCell2/dev_data/sim_list_noise_sd1090523.rds")

sim_list_noise_sd1 <- makeSimulations(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = TRUE, noise_sd = 1)
# saveRDS(sim_list_noise_sd1, "/bigdata/almogangel/xCell2/dev_data/sim_list_noise_sd1090523.rds")

sim_list_noise_sd1.5 <- makeSimulations(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = TRUE, noise_sd = 1.5)
# saveRDS(sim_list_noise_sd1.5, "/bigdata/almogangel/xCell2/dev_data/sim_list_noise_sd1090523.rds")

sim_list_noise_sd2 <- makeSimulations(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = TRUE, noise_sd = 2)
# saveRDS(sim_list_noise_sd1, "/bigdata/almogangel/xCell2/dev_data/sim_list_noise_sd1090523.rds")

sim_list_noise_adjusted <- makeSimulations(ref, labels, mixture_fractions, dep_list, n_ct_sim = 20, add_noise = TRUE)
# saveRDS(sim_list_noise_adjusted, "/bigdata/almogangel/xCell2/dev_data/sim_list_noise_sd1090523.rds")




# Investigate best parameters of probs, diff and n_genes for each cell type using the simulation
# Compare the results with the validation datasets

getSimulationCors <- function(sim, ctoi, signatures_collection){

  fractions <- as.numeric(colnames(sim))
  signatures_ctoi <- signatures_collection[startsWith(names(signatures_collection), paste0(ctoi, "#"))]
  sim_ranked <- singscore::rankGenes(sim)

  scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
    singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })
  scores <- t(scores)
  rownames(scores) <- names(signatures_ctoi)

  sig_cors <- apply(scores, 1, function(sig){
    cor(sig, fractions, method = "spearman")
  })

  return(sig_cors)

}

makeSimulationResultsTidy <- function(simulations_results, simulation_type){
  enframe(simulations_results, name = "celltype", value = "simulations") %>%
    unnest_longer(simulations, indices_to = "simulation_rep") %>%
    unnest_longer(simulations, indices_to = "signature", values_to = "cor") %>%
    group_by(celltype, signature) %>%
    summarise(median_cor = median(cor)) %>%
    separate(signature, sep = "#_", into = c("remove", "signature_info"), remove = FALSE) %>%
    select(-remove) %>%
    separate(signature_info, sep = "_", into = c("probs", "diff", "n_genes"), remove = TRUE) %>%
    mutate(test = simulation_type) %>%
    ungroup() %>%
    return(.)
}

# Get simulations correlation scores
celltypes <- names(sim_list)
simulations_results <- lapply(celltypes, function(ctoi){
  print(ctoi)

  pbapply::pblapply(sim_list[[ctoi]], function(sim){
    getSimulationCors(sim, ctoi, signatures_collection)
  })

})
names(simulations_results) <- celltypes
# saveRDS(simulations_results, "/bigdata/almogangel/xCell2/dev_data/simulations_results080523.rds")

celltypes <- names(sim_list_noise_sd0.5)
simulations_results_noise_sd0.5 <- lapply(celltypes, function(ctoi){
  print(ctoi)

  pbapply::pblapply(sim_list_noise_sd0.5[[ctoi]], function(sim){
    getSimulationCors(sim, ctoi, signatures_collection)
  })

})
names(simulations_results_noise_sd0.5) <- celltypes
# saveRDS(simulations_results_noise_sd0.5, "/bigdata/almogangel/xCell2/dev_data/simulations_results_noise_sd1080523.rds")

celltypes <- names(sim_list_noise_sd1)
simulations_results_noise_sd1 <- lapply(celltypes, function(ctoi){
  print(ctoi)

  pbapply::pblapply(sim_list_noise_sd1[[ctoi]], function(sim){
    getSimulationCors(sim, ctoi, signatures_collection)
  })

})
names(simulations_results_noise_sd1) <- celltypes
# saveRDS(simulations_results_noise_sd1, "/bigdata/almogangel/xCell2/dev_data/simulations_results_noise_sd1080523.rds")

celltypes <- names(sim_list_noise_sd1.5)
simulations_results_noise_sd1.5 <- lapply(celltypes, function(ctoi){
  print(ctoi)

  pbapply::pblapply(sim_list_noise_sd1.5[[ctoi]], function(sim){
    getSimulationCors(sim, ctoi, signatures_collection)
  })

})
names(simulations_results_noise_sd1.5) <- celltypes
# saveRDS(simulations_results_noise_sd1.5, "/bigdata/almogangel/xCell2/dev_data/simulations_results_noise_sd1080523.rds")

celltypes <- names(sim_list_noise_sd2)
simulations_results_noise_sd2 <- lapply(celltypes, function(ctoi){
  print(ctoi)

  pbapply::pblapply(sim_list_noise_sd2[[ctoi]], function(sim){
    getSimulationCors(sim, ctoi, signatures_collection)
  })

})
names(simulations_results_noise_sd2) <- celltypes
# saveRDS(simulations_results_noise_sd2, "/bigdata/almogangel/xCell2/dev_data/simulations_results_noise_sd1080523.rds")

celltypes <- names(sim_list_noise_adjusted)
simulations_results_noise_adjusted <- lapply(celltypes, function(ctoi){
  print(ctoi)

  pbapply::pblapply(sim_list_noise_adjusted[[ctoi]], function(sim){
    getSimulationCors(sim, ctoi, signatures_collection)
  })

})
names(simulations_results_noise_adjusted) <- celltypes
# saveRDS(simulations_results_noise_adjusted, "/bigdata/almogangel/xCell2/dev_data/simulations_results_noise_sd1080523.rds")


# Analyze simulations
simulations_results_tidy <- makeSimulationResultsTidy(simulations_results, simulation_type = "simulation - no noise")
simulations_results_noise_sd0.5_tidy <- makeSimulationResultsTidy(simulations_results_noise_sd0.5, simulation_type = "simulation - 0.5sd noise")
simulations_results_noise_sd1_tidy <- makeSimulationResultsTidy(simulations_results_noise_sd1, simulation_type = "simulation - 1sd noise")
simulations_results_noise_sd1.5_tidy <- makeSimulationResultsTidy(simulations_results_noise_sd1.5, simulation_type = "simulation - 1.5sd noise")
simulations_results_noise_sd2_tidy <- makeSimulationResultsTidy(simulations_results_noise_sd2, simulation_type = "simulation - 2sd noise")
simulations_results_noise_adjusted_tidy <- makeSimulationResultsTidy(simulations_results_noise_adjusted, simulation_type = "simulation - adjusted noise")


validation_results <- tibble(id = names(all_ds_cors), cor = all_ds_cors)

validation_results_tidy <- validation_results %>%
  separate(id, sep = "#", into = c("dataset", "signature"), extra = "merge") %>%
  separate(signature, sep = "#", into = "celltype", remove = FALSE, extra = "drop") %>%
  group_by(celltype, signature) %>%
  summarise(median_cor = median(cor)) %>%
  separate(signature, sep = "#_", into = c("remove", "signature_info"), remove = FALSE) %>%
  select(-remove) %>%
  separate(signature_info, sep = "_", into = c("probs", "diff", "n_genes"), remove = TRUE) %>%
  mutate(test = "validation") %>%
  ungroup()

top_sigs_no_noise <- simulations_results_tidy %>%
  group_by(celltype) %>%
  top_frac(n=0.1, wt=median_cor) %>%
  pull(signature)

top_sigs_noise_0.5sd <- simulations_results_noise_sd0.5_tidy %>%
  group_by(celltype) %>%
  top_frac(n=0.1, wt=median_cor) %>%
  pull(signature)

top_sigs_noise_1sd <- simulations_results_noise_sd1_tidy %>%
  group_by(celltype) %>%
  top_frac(n=0.1, wt=median_cor) %>%
  pull(signature)

top_sigs_noise_1.5sd <- simulations_results_noise_sd1.5_tidy %>%
  group_by(celltype) %>%
  top_frac(n=0.1, wt=median_cor) %>%
  pull(signature)

top_sigs_noise_2sd <- simulations_results_noise_sd2_tidy %>%
  group_by(celltype) %>%
  top_frac(n=0.1, wt=median_cor) %>%
  pull(signature)

top_sigs_noise_adjusted <- simulations_results_noise_adjusted_tidy %>%
  group_by(celltype) %>%
  top_frac(n=0.1, wt=median_cor) %>%
  pull(signature)


validation_results_tidy_top_sigs_no_noise <- validation_results_tidy %>%
  filter(signature %in% top_sigs_no_noise) %>%
  mutate(test = "validation - top 10% sigs in simulation - no noise")

validation_results_tidy_top_sigs_noise_0.5sd <- validation_results_tidy %>%
  filter(signature %in% top_sigs_noise_0.5sd) %>%
  mutate(test = "validation - top 10% sigs in simulation - 0.5sd noise")

validation_results_tidy_top_sigs_noise_1sd <- validation_results_tidy %>%
  filter(signature %in% top_sigs_noise_1sd) %>%
  mutate(test = "validation - top 10% sigs in simulation - 1sd noise")

validation_results_tidy_top_sigs_noise_1.5sd <- validation_results_tidy %>%
  filter(signature %in% top_sigs_noise_1.5sd) %>%
  mutate(test = "validation - top 10% sigs in simulation - 1.5sd noise")

validation_results_tidy_top_sigs_noise_2sd <- validation_results_tidy %>%
  filter(signature %in% top_sigs_noise_2sd) %>%
  mutate(test = "validation - top 10% sigs in simulation - 2sd noise")

validation_results_tidy_top_sigs_noise_adjusted <- validation_results_tidy %>%
  filter(signature %in% top_sigs_noise_adjusted) %>%
  mutate(test = "validation - top 10% sigs in simulation")

shared_celltypes <- intersect(unique(simulations_results_tidy$celltype), unique(validation_results_tidy$celltype))
celltypes2use <- c("T cell", "monocyte", "CD4-positive, alpha-beta T cell", "memory T cell", "B cell", "eosinophil", "natural killer cell", "CD8-positive, alpha-beta T cell",
                   "neutrophil", "helper T cell")

rbind(validation_results_tidy, validation_results_tidy_top_sigs_noise_adjusted) %>%
  filter(celltype %in% celltypes2use) %>%
  ggplot(., aes(x=test, y=median_cor, fill=test)) +
  geom_boxplot() +
  facet_wrap(~celltype, scales = "free_x")
