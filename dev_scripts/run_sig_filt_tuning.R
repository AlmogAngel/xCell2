library(tidyverse)
library(xCell2)



setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")


# "/bigdata/almogangel/xCell2/dev_scripts/prep_ref_val_pairs.R"
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
# sc.refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val.rds")
# Load validation data
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
# sc.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc.vals.rds")

refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))


# Load references matrices
refsRDSList <- lapply(refList, function(ref_type){
  refs <- lapply(ref_type, function(ref){
    # Load reference
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
    ref.in
  })
  names(refs) <- ref_type
  refs
})

vals.refs.res <- refval.tbl %>%
  # Get number of samples in the validation dataset
  mutate(n_val_samples = ncol(cyto.vals$truth[[val_type]][[val_dataset]])) %>%
  filter(n_shared_celltypes > 2) %>%
  mutate(method = "xCell2", .before = everything())



# xCell 2.0 settings
thisseed <- 123

# Signatures generation tuning
# TODO: DO NOT FORGET TO RUN use_sim2filter (TRUE/FALSE)
# max_rho_cutoff = c(0.3, 0.5)
# top_frac_sigs_ds = c(0.25, 0.5, 0.75)
# min_ds_frac = c(0.25, 0.5, 0.75)
# top_sigs_frac = c(0.25, 0.5, 0.75)
#
# grid <- expand.grid(max_rho_cutoff = max_rho_cutoff, top_frac_sigs_ds = top_frac_sigs_ds, min_ds_frac = min_ds_frac, top_sigs_frac = top_sigs_frac)

# write.table(grid, "/bigdata/almogangel/xCell2_data/benchmarking_data/sig_filt_tuning_1/grid_table_23jun.txt", sep = "\t")
# grid 8:
# max_rho_cutoff = 0.5
# top_frac_sigs_ds = 0.25
# min_ds_frac = 0.5
# top_sigs_frac = 0.25



min_top_genes_frac = c(0.25, 0.5, 0.75)
essen_gene_cutoff = c(0.1, 0.25, 0.5, 0.75)
min_filt_sigs = c(3, 5, 10)
top_sigs_frac = c(0.05, 0.1, 0.25)

grid <- expand.grid(min_top_genes_frac = min_top_genes_frac, essen_gene_cutoff = essen_gene_cutoff,
                    min_filt_sigs = min_filt_sigs, top_sigs_frac = top_sigs_frac)

# write.table(grid, "/bigdata/almogangel/xCell2_data/benchmarking_data/sig_filt_tuning_2/grid_table_25jun.txt", sep = "\t")


for (j in 1:nrow(grid)) {

  file <-  paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/sig_filt_tuning_2/grid_", j, "_res.rds")


  if (file.exists(file)) {
    print(paste0("Gird: ", j, " exist..."))
    next
  }else{
    print(paste0("Gird: ", j, "/", nrow(grid)))
  }

  xCell2results <- parallel::mclapply(1:nrow(vals.refs.res), function(i){

    print(paste0("-------------------- ", i, "/", nrow(vals.refs.res), " --------------------"))

    # Load data
    val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
    print(val_ref)
    mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
    ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
    labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
    lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
    refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)
    valType <- vals.refs.res[i,]$val_type
    valDataset <- vals.refs.res[i,]$val_dataset
    refName <- vals.refs.res[i,]$ref_name
    valDataType <- vals.refs.res[i,]$val_data_type[[1]]


    # Load filtering data
    if (valDataType == "array") {
      filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/array_filtering_data.rds")
    }else{
      filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/rnaseq_filtering_data.rds")
    }


    # Remove current validation from filtering data
    filt_datasets <- gsub(x = names(filtering_data$mixture), pattern = "#.*", replacement = "")
    filtering_data$mixture <- filtering_data$mixture[filt_datasets != valDataset]
    filtering_data$truth <- filtering_data$truth[filt_datasets != valDataset]

    xcell2_object <- xCell2::xCell2Train(ref = ref.in, labels = labels, mix = mix.in, filtering_data = filtering_data, ref_type = refType, lineage_file = lineage_file,
                                         num_threads = 25, return_analysis = FALSE, use_sillover = FALSE, use_ontology = TRUE, return_signatures = FALSE,
                                         max_rho_cutoff = 0.5, top_frac_sigs_ds = 0.25,
                                         min_ds_frac = 0.5,
                                         use_sim2filter = FALSE,
                                         min_top_genes_frac = grid$min_top_genes_frac[j],
                                         essen_gene_cutoff = grid$essen_gene_cutoff[j],
                                         min_filt_sigs = grid$min_filt_sigs[j],
                                         top_sigs_frac = grid$top_sigs_frac[j])
                                        # TODO: tune use_sim2filter -> TRUE




    res <- xCell2::xCell2Analysis(mix.in, xcell2object = xcell2_object, raw_scores = TRUE,
                                  spillover = FALSE, spillover_alpha = 0, num_threads = 1)

    return(res)


  }, mc.cores = 10)

  saveRDS(xCell2results, file)

}

# Analyze grid

if (FALSE) {

  vals.refs <- vals.refs.res %>% select(val_dataset, ref_name, n_val_samples, ref_type, val_data_type)

  res.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/sig_gen_tuning_2/grid2_1_res.rds")
  vals.refs$res <- res.in

  res.in <- vals.refs %>%
    unnest(res)

  res.in <- res.in %>%
    unnest_longer(sigs_cors, values_to = "sig_cor", indices_to = "sig_name") %>%
    separate(sig_name, into = c("ct", "sig_prob", "sig_diff", "n_sig_genes", "n_ct_passed", "frac_ct_passed"), sep = "_") %>%
    select(-c(ct, mean_cor, n_ct_passed))


  res.in <- res.in %>%
    mutate(sig_prob = as.numeric(sig_prob),
           sig_diff = as.numeric(sig_diff),
           n_sig_genes = as.numeric(n_sig_genes),
           frac_ct_passed = as.numeric(frac_ct_passed))


  res.in <- res.in %>%
    mutate(sig_cor = ifelse(sig_cor == 1, 0.99999, sig_cor)) %>%
    mutate(z = 0.5 * log((1 + sig_cor) / (1 - sig_cor))) %>%
    mutate(z_weighted = log(n_val_samples) * z)

  data <- res.in

  data$val_dataset <- as.factor(data$val_dataset)
  data$ref_name <- as.factor(data$ref_name)
  data$celltype <- as.factor(data$celltype)
  data$ref_type <- as.factor(data$ref_type)
  data$val_data_type <- as.factor(data$val_data_type)

  # Create bins for each parameter
  data <- data %>%
    mutate(
      n_sig_genes_bin = cut(n_sig_genes, breaks = seq(0, max(n_sig_genes, na.rm = TRUE), by = 10), include.lowest = TRUE),
      sig_prob = factor(sig_prob),
      sig_diff = factor(sig_diff),
      ref_type = factor(ref_type),
      val_data_type = factor(val_data_type),
      frac_ct_passed_bin = cut(frac_ct_passed, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
    )

  model <- lme4::lmer(z_weighted ~ sig_prob + sig_diff + n_sig_genes_bin + frac_ct_passed_bin + (1 | val_dataset)  + (1 | ref_name)  + (1 | celltype),  data = data)


  # Marginal effects for sig_prob
  effects_sig_prob <- ggeffects::ggpredict(model, terms = "sig_prob")
  plot(effects_sig_prob) + labs(title = "Marginal Effects of sig_prob", x = "sig_prob", y = "Predicted z_weighted")

  # Marginal effects for sig_diff
  effects_sig_diff <- ggeffects::ggpredict(model, terms = "sig_diff")
  plot(effects_sig_diff) + labs(title = "Marginal Effects of sig_diff", x = "sig_diff", y = "Predicted z_weighted")

  # Marginal effects for n_sig_genes
  effects_n_sig_genes <- ggeffects::ggpredict(model, terms = "n_sig_genes_bin")
  plot(effects_n_sig_genes) + labs(title = "Marginal Effects of n_sig_genes_bin", x = "n_sig_genes_bin", y = "Predicted z_weighted")

  # Marginal effects for frac_ct_passed
  effects_frac_ct_passed <- ggeffects::ggpredict(model, terms = "frac_ct_passed_bin")
  plot(effects_frac_ct_passed) + labs(title = "Marginal Effects of frac_ct_passed_bin", x = "frac_ct_passed_bin", y = "Predicted z_weighted")


  model_interaction_prob_diff <- lme4::lmer(z_weighted ~ sig_prob * sig_diff + (1 | val_dataset) + (1 | ref_name) + (1 | celltype), data = data)

  summary(model_interaction_prob_diff)

  effects_sig_prob_sig_diff <- ggeffects::ggpredict(model_interaction_prob_diff, terms = c("sig_prob", "sig_diff"))
  plot(effects_sig_prob_sig_diff) + labs(title = "Interaction: sig_prob and sig_diff", x = "sig_prob", y = "Predicted z_weighted")

  # Convert the effects to a data frame for ggplot2
  effects_df <- as.data.frame(effects_sig_prob_sig_diff)
  colnames(effects_df)[1] <- "sig_prob"
  colnames(effects_df)[6] <- "sig_diff"

  color_palette <- colorRampPalette(c("darkred" ,"red", "white", "lightblue1", "darkblue"))(50)


  ggplot(effects_df, aes(x = sig_prob, y = sig_diff, fill = predicted)) +
    geom_tile() +
    scale_fill_gradientn(colors = color_palette) +
    labs(title = "Interaction Effect of sig_prob and sig_diff on z_weighted",
         x = "sig_prob",
         y = "sig_diff",
         fill = "Predicted z_weighted") +
    theme_minimal()


  model_interaction_genes_frac <- lme4::lmer(z_weighted ~ n_sig_genes_bin * frac_ct_passed_bin + (1 | val_dataset) + (1 | ref_name) + (1 | celltype), data = data)

  # Summarize the model
  summary(model_interaction_genes_frac)

  # Get the predicted values for the interaction between n_sig_genes_bin and frac_ct_passed_bin
  effects_genes_frac <- ggpredict(model_interaction_genes_frac, terms = c("n_sig_genes_bin", "frac_ct_passed_bin"))

  # Convert the effects to a data frame for ggplot2
  effects_df_genes_frac <- as.data.frame(effects_genes_frac)
  colnames(effects_df_genes_frac)[1] <- "n_sig_genes_bin"
  colnames(effects_df_genes_frac)[6] <- "frac_ct_passed_bin"



  # Define the custom color palette
  color_palette <- colorRampPalette(c("darkred", "red", "white", "lightblue1", "darkblue"))(50)

  # Create a heatmap with custom color palette
  ggplot(effects_df_genes_frac, aes(x = n_sig_genes_bin, y = frac_ct_passed_bin, fill = predicted)) +
    geom_tile() +
    scale_fill_gradientn(colors = color_palette) +
    labs(title = "Interaction Effect of n_sig_genes_bin and frac_ct_passed_bin on z_weighted",
         x = "n_sig_genes_bin",
         y = "frac_ct_passed_bin",
         fill = "Predicted z_weighted") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))



  # Fit the model including the interaction between sig_prob, sig_diff, ref_type, and val_data_type
  model_interaction_prob_diff_data <- lme4::lmer(z_weighted ~ sig_prob * sig_diff * ref_type * val_data_type + (1 | val_dataset) + (1 | ref_name) + (1 | celltype), data = data)

  # Summarize the model
  summary(model_interaction_prob_diff_data)

  # Get the predicted values for the interaction between sig_prob, sig_diff, ref_type, and val_data_type
  effects_prob_diff_data <- ggpredict(model_interaction_prob_diff_data, terms = c("sig_prob", "sig_diff", "ref_type", "val_data_type"))

  # Convert the effects to a data frame for ggplot2
  effects_df_prob_diff_data <- as.data.frame(effects_prob_diff_data)
  colnames(effects_df_prob_diff_data)[1] <- "sig_prob"
  colnames(effects_df_prob_diff_data)[6] <- "sig_diff"
  colnames(effects_df_prob_diff_data)[7] <- "ref_type"
  colnames(effects_df_prob_diff_data)[8] <- "val_data_type"

  # Define the custom color palette
  color_palette <- colorRampPalette(c("darkred", "red", "white", "lightblue1", "darkblue"))(50)

  # Create a heatmap with custom color palette for each ref_type and val_data_type
  ggplot(effects_df_prob_diff_data, aes(x = sig_prob, y = sig_diff, fill = predicted)) +
    geom_tile() +
    scale_fill_gradientn(colors = color_palette) +
    facet_grid(ref_type ~ val_data_type) +
    labs(title = "Interaction Effect of sig_prob and sig_diff on z_weighted by Data Type",
         x = "sig_prob",
         y = "sig_diff",
         fill = "Predicted z_weighted") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

}







