validateInputs <- function(ref, labels, data_type){
  if (length(unique(labels$label)) < 3) {
    stop("Reference must have at least 3 cell types")
  }

  if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
    stop("ref must be one of those classes: matrix, dgCMatrix, Matrix")
  }

  if (!"data.frame" %in% class(labels)) {
    stop("labels must be a dataframe.")
  }

  if (!data_type %in% c("rnaseq", "array", "sc")) {
    stop("data_type should be 'rnaseq', 'array' or 'sc'.")
  }

  if (sum(grepl("_", labels$label)) != 0) {
    message("Changing underscores to dashes in cell-types labels!")
    labels$label <- gsub("_", "-", labels$label)
  }


  out <- list(ref = ref,
              labels = labels)
  return(out)

}
prepareRef <- function(ref, data_type){

  if (data_type == "sc") {
    if(max(ref) >= 50){
      message("Normalizing and transforming scRNA-Seq reference to log1p-space (maximum expression value >= 50).")
      genes_names <- rownames(ref)
      ref.srt <- Seurat::CreateSeuratObject(counts = ref)
      ref.srt <- Seurat::NormalizeData(ref.srt, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
      ref.norm <-  ref.srt@assays$RNA@data

      # Check if Seurat changed genes names
      if (!all(rownames(ref.norm) %in% genes_names)) {
        # rownames(ref.norm) <- genes_names # Because Seurat change genes names from "_" to "-"
        stop("Seurat genes name error")
      }

      return(ref.norm)
    }else{
      message("Assuming reference is already in log1p-space (maximum expression value < 50).")
      return(ref)
    }

  }else{

    if(max(ref) >= 50){
      message("Transforming reference to log2-space (maximum expression value >= 50).")
      ref.norm <- log2(ref+1)
      return(ref.norm)
    }else{
      message("Assuming reference is already in log2-space (maximum expression value < 50).")
      return(ref)
    }

  }

}
makePureCTMat <- function(ref, labels, use_median){

  celltypes <- unique(labels$label)

  pure_ct_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      if(use_median){
        type_vec <- if("matrix" %in% class(ref)) Rfast::rowMedians(ref[,type_samples]) else sparseMatrixStats::rowMedians(ref[,type_samples])
      }else{
        type_vec <- if("matrix" %in% class(ref)) Rfast::rowmeans(ref[,type_samples]) else Matrix::rowMeans(ref[,type_samples])
      }
    }
  })
  rownames(pure_ct_mat) <- rownames(ref)

  return(pure_ct_mat)
}
getCellTypeCorrelation <- function(pure_ct_mat, data_type){

  celltypes <- colnames(pure_ct_mat)

  if (data_type != "sc") {

    # Use top 10% most variable genes
    genes_var <- apply(pure_ct_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.9, na.rm=TRUE)
    pure_ct_mat <- pure_ct_mat[genes_var > most_var_genes_cutoff,]

  }else{

    # Use top 1% most variable genes
    genes_var <- apply(pure_ct_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.99, na.rm=TRUE)
    pure_ct_mat <- pure_ct_mat[genes_var > most_var_genes_cutoff,]

  }


  # Make correlation matrix
  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

  # TODO: Change for loop to apply function to measure time
  for (i in 1:nrow(lower_tri_coord)) {
    celltype_i <- rownames(cor_mat)[lower_tri_coord[i, 1]]
    celltype_j <- colnames(cor_mat)[lower_tri_coord[i, 2]]
    cor_mat[lower_tri_coord[i, 1], lower_tri_coord[i, 2]] <- cor(pure_ct_mat[,celltype_i], pure_ct_mat[,celltype_j], method = "spearman")
    cor_mat[lower_tri_coord[i, 2], lower_tri_coord[i, 1]] <- cor(pure_ct_mat[,celltype_i], pure_ct_mat[,celltype_j], method = "spearman")
  }

  return(cor_mat)
}
getDependencies <- function(lineage_file_checked){
  ont <- read_tsv(lineage_file_checked, show_col_types = FALSE) %>%
    mutate_all(as.character)

  celltypes <- pull(ont[,2])
  celltypes <- gsub("_", "-", celltypes)
  dep_list <- vector(mode = "list", length = length(celltypes))
  names(dep_list) <- celltypes

  for (i in 1:nrow(ont)) {
    descendants <-  gsub("_", "-", strsplit(pull(ont[i,3]), ";")[[1]])
    descendants <- descendants[!is.na(descendants)]

    ancestors <-  gsub("_", "-", strsplit(pull(ont[i,4]), ";")[[1]])
    ancestors <- ancestors[!is.na(ancestors)]

    dep_list[[i]] <- list("descendants" = descendants, "ancestors" = ancestors)

  }

  return(dep_list)
}
makeQuantiles <- function(ref, labels, probs, dep_list, include_descendants){

  celltypes <- unique(labels[,2])

  quantiles_mat_list <-  pbapply::pblapply(celltypes, function(type){

    if (include_descendants) {
      # Include all the descendants of the cell type for quantiles calculations
      descen_cells <- dep_list[[type]]$descendants
      type_samples <- labels[,2] == type | labels[,2] %in% descen_cells
    }else{
      type_samples <- labels[,2] == type
    }


    if (sum(type_samples) == 1) {
      # If there is one sample for this cell type -> duplicate the sample to make a data frame
      type.df <- cbind(ref[,type_samples], ref[,type_samples])
    }else{
      type.df <- ref[,type_samples]
    }

    # Calculate quantiles
    # TODO: Balance quantiles by dataset
    type_quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1-probs))), na.rm=TRUE))
  })
  names(quantiles_mat_list) <- celltypes

  return(quantiles_mat_list)
}
createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Set signature parameters
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)

    # Generate signatures
    type_sigs <- list()
    for(i in 1:nrow(param.df)){

      # Get a Boolean matrices with genes that pass the quantiles criteria
      diff <- param.df[i, ]$diff_vals # difference criteria
      lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile criteria
      upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile criteria
      diff_genes.mat <- sapply(not_dep_celltypes, function(x){
        get(type, quantiles_matrix)[lower_prob,] > get(x, quantiles_matrix)[upper_prob,] + diff
      })


      if(weight_genes){
        # Score genes using cell types correlations as weights
        type_weights <- cor_mat[type, not_dep_celltypes]
        type_weights[type_weights < 0.001] <- 0.001 # Fix minimum correlation to avoid zero and negative correlations
        gene_scores <- apply(diff_genes.mat, 1, function(x){
          sum(type_weights[which(x)])
        })

      }else{
        # All cell types scores are the same
        gene_scores <- apply(diff_genes.mat, 1, function(x){
          sum(x)
        })

      }


      gene_passed <- sort(gene_scores[gene_scores > 0], decreasing = TRUE)

      # If less than min_genes passed move to next parameters
      if (length(gene_passed) < min_genes) {
        next
      }


      # Round and sort genes scores
      top_scores <- sort(unique(round(gene_passed-0.5)), decreasing = TRUE)

      # Take top 10 highest scores
      # TODO: test different n_tops parameters (currently 10 is default)
      n_top <- ifelse(length(top_scores) > 10, 10, length(top_scores))

      for (score in top_scores[1:n_top]) {

        n_genes <- sum(gene_passed >= score)

        if (n_genes < min_genes) {
          next
        }

        if (n_genes > max_genes) {
          break
        }

        sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, sep = "_")
        type_sigs[[sig_name]] <- names(which(gene_passed >= score))
      }

    }

    if (length(type_sigs) == 0) {
      warning(paste0("No signatures found for ", type))
    }


    # Remove duplicate signatures
    type_sigs_sorted <- lapply(type_sigs, function(x) sort(x))
    type_sigs_sorted_collapsed <- sapply(type_sigs_sorted, paste, collapse = ",")
    duplicated_sigs <- duplicated(type_sigs_sorted_collapsed)
    type_sigs <- type_sigs[!duplicated_sigs]


    return(type_sigs)
  }


  celltypes <- unique(labels[,2])

  all_sigs <- pbapply::pblapply(celltypes, function(type){
    getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes)
  })
  all_sigs <- unlist(all_sigs, recursive = FALSE)


  if (length(all_sigs) == 0) {
    warning("No signatures found for reference!")
  }


  return(all_sigs)
}
makeSimulations <- function(ref, labels, pure_ct_mat, cor_mat, dep_list, sim_fracs, n_sims, n_samples_sim, add_noise, seed = 123, simple){

  set.seed(seed)
  celltypes <- unique(labels$label)

  makeSamplesPool <- function(labels, ctoi){

    # Mix CTOI samples by datasets
    ctoi_samples_pool <- c()
    while(!all(labels[labels$label == ctoi,]$sample %in% ctoi_samples_pool)) {
      ctoi_samples_pool <- c(ctoi_samples_pool,
                             labels %>%
                               filter(label == ctoi & !sample %in% ctoi_samples_pool) %>%
                               slice_head(n = 1, by = dataset) %>%
                               pull(sample))
    }

    return(ctoi_samples_pool)
  }
  makeFractionMatrixCTOI <- function(ref, sim_fracs, ctoi_samples2use){

    if (length(ctoi_samples2use) == 1) {
      ctoi_mean_expression <- ref[,ctoi_samples2use]
    }else{
      ctoi_mean_expression <- if("matrix" %in% class(ref)) Rfast::rowmeans(ref[,ctoi_samples2use]) else Matrix::rowMeans(ref[,ctoi_samples2use])
    }

    ctoi_frac_mat <- matrix(rep(ctoi_mean_expression, length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs)) %*% diag(sim_fracs)
    rownames(ctoi_frac_mat) <- rownames(ref)

    return(ctoi_frac_mat)
  }
  makeFractionMatrixControls <- function(ref, labels, sim_fracs, dep_list, controls, n_samples_sim, max_cts_controls){

    max_cts_controls <-  ifelse(max_cts_controls > length(controls), length(controls), max_cts_controls)

    if (length(controls) > 1) {
      # Sample independent controls
      controls_shufled <- sample(controls)
      controls2use <- c()
      deps_seen <- c()
      for (ctrl in controls_shufled) {

        if (length(controls2use) == max_cts_controls) {
          break
        }

        deps_seen <- c(deps_seen, unname(unlist(dep_list[[ctrl]])))

        if (ctrl %in% deps_seen) {
          next
        }

        controls2use <- c(ctrl, controls2use)

      }
    }else{
      controls2use <- controls
    }


    # Generate controls matrix
    controls_mat <- sapply(controls2use, function(ctrl){

      control_samples <- labels[labels$label == ctrl,]$sample

      if (is.null(n_samples_sim)) {
        # For simple simulation take all control samples
        samples2use <- control_samples
      }else{
        # For complex simulation pick random control samples by dataset
        samples2use <- c()
        if (length(control_samples) <= n_samples_sim) {
          samples2use <- control_samples
        }else{
          while(length(samples2use) < n_samples_sim){
            samples2use <-labels %>%
              filter(label == ctrl & !sample %in% samples2use) %>%
              group_by(dataset) %>%
              slice_sample(n=1) %>%
              pull(sample) %>%
              c(samples2use, .)
          }
        }
      }


      if (length(samples2use) == 1) {
        ref[,samples2use]
      }else{
        if("matrix" %in% class(ref)) Rfast::rowmeans(ref[,samples2use]) else Matrix::rowMeans(ref[,samples2use])
      }

    })
    rownames(controls_mat) <- rownames(ref)


    # Multiple control mixture by fractions
    generateFractions <- function(target_sum, n_fracs = length(controls2use)) {
      # Generate random numbers from a uniform distribution
      numbers <- runif(n_fracs)

      # Scale the numbers so that their sum equals target_sum
      fracs <- numbers / sum(numbers) * target_sum

      return(fracs)
    }

    controls_frac_mat <- sapply(1-sim_fracs, function(frac){
      fracs <- generateFractions(target_sum = frac)
      if (length(fracs) == 1) {
        controls_mat * fracs
      }else{
        if("matrix" %in% class(ref)) Rfast::rowmeans(controls_mat %*% diag(fracs)) else Matrix::rowMeans(controls_mat %*% diag(fracs))
      }
    })
    rownames(controls_frac_mat) <- rownames(ref)

    return(controls_frac_mat)
  }


  sim_list <- pbapply::pblapply(celltypes, function(ctoi){

    # Generate CTOI samples pool
    ctoi_samples_pool <- makeSamplesPool(labels, ctoi)

    # Number of CTOI samples to use
    if (!simple) {
      n_samples_sim <- ifelse(n_samples_sim < length(ctoi_samples_pool), n_samples_sim, length(ctoi_samples_pool))
    }

    # Get control(s) cell types
    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
    controls <- celltypes[!celltypes %in% dep_cts]
    if (simple) {
      # Pick one control for the simple simulation
      controls <- names(sort(cor_mat[ctoi, controls])[1])
    }

    # Get number of available datasets of this cell type (for adding noise)
    n_ctoi_ds <- length(unique(labels[labels$label == ctoi,]$dataset))

    # Generate n_sims simulations
    ctoi_sim_list <- lapply(1:n_sims, function(i){

      # Make CTOI fraction matrix
      if (is.null(n_samples_sim)) {
        samples2use <- ctoi_samples_pool
      }else{
        samples2use <- ctoi_samples_pool[1:n_samples_sim]
      }
      ctoi_frac_mat <- makeFractionMatrixCTOI(ref, sim_fracs, ctoi_samples2use = samples2use)

      # Move samples2use to the end of the vector
      ctoi_samples_pool <- c(ctoi_samples_pool[!ctoi_samples_pool %in% samples2use], samples2use)

      # Make control(s) fractions matrix
      if (simple) {
        controls_frac_mat <- makeFractionMatrixControls(ref, labels, sim_fracs, dep_list, controls, n_samples_sim, max_cts_controls = 1)
        colnames(controls_frac_mat) <- paste0(controls, "%%", 1-sim_fracs)
      }else{
        controls_frac_mat <- makeFractionMatrixControls(ref, labels, sim_fracs, dep_list, controls, n_samples_sim, max_cts_controls = 8)
      }

      # Combine CTOI and control(s) fractions matrix
      simulation <- ctoi_frac_mat + controls_frac_mat
      colnames(simulation) <- paste0(ctoi, "%%", sim_fracs)

      # Add noise
      if (add_noise) {
        noise_sd <- 1/n_ctoi_ds
        # noise_sd <- 1
        noise <- matrix(rnorm(nrow(simulation) * ncol(simulation), mean = 0, sd = noise_sd),
                        nrow = nrow(simulation), ncol = ncol(simulation))
        simulation <- simulation + noise
        simulation <- pmax(simulation, 0)
      }

      if (simple) {
        # Add another sample to estimate the control background expression for the simple simulation
        simulation <- cbind(simulation, controls_frac_mat[,1])
        colnames(simulation)[ncol(simulation)] <- colnames(controls_frac_mat)[1]
      }

      simulation

    })

    if (n_sims == 1) {
      ctoi_sim_list[[1]]
    }else{
      names(ctoi_sim_list) <- paste0("sim-", 1:n_sims)
      ctoi_sim_list
    }
  })
  names(sim_list) <- celltypes

  return(sim_list)

}
filterSignatures <- function(sim, signatures_collection, dep_list, n_null = 1000, n_cpu = 10){

  scoreSim <- function(sim, signatures_collection, dep_list, n_null){

    # Rank mixture
    sim_ranked <- singscore::rankGenes(sim)

    # Generate null distribution for p-values
    sigs_genes <- unique(unlist(signatures_collection))
    non_sigs_genes <- rownames(sim_ranked)[!rownames(sim_ranked) %in% sigs_genes]
    sigs_lengths <- unique(lengths(signatures_collection))
    all_lengths_null_scores <- pbapply::pblapply(sigs_lengths, function(len){
      tmp_genes <- sample(non_sigs_genes, len)
      singscore::generateNull(
        upSet = tmp_genes,
        rankData = sim_ranked,
        subSamples = 1:ncol(sim_ranked),
        centerScore = FALSE,
        B = n_null,
        ncores = n_cpu,
        seed = 1)
    })
    names(all_lengths_null_scores) <- sigs_lengths


    # Score signatures and get p-values
    sigs_scores_pvals <- pbapply::pblapply(signatures_collection, function(sig){
      sig_len <- length(sig)
      scoredf <- singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)
      pvals <- singscore::getPvals(all_lengths_null_scores[[as.character(sig_len)]], scoredf)
      tibble("sim_name" = names(pvals), "score" = scoredf$TotalScore, "pval" = pvals)
    })
    names(sigs_scores_pvals) <- names(signatures_collection)


    # Make results tidy
    sigs_scores_pvals_tidy <- enframe(sigs_scores_pvals, name = "signature") %>%
      unnest(value) %>%
      separate(sim_name, into = c("sim_celltype", "sim_control"), sep = "%%") %>%
      separate(signature, into = "sig_celltype", sep = "#", extra = "drop", remove = FALSE)

    # Clean scores
    sigs_scores_pvals_tidy_clean <- sigs_scores_pvals_tidy %>%
      filter(sig_celltype != sim_control) %>% # Signature cell type cannot be the same as the control
      rowwise() %>%
      filter(!sim_celltype %in% unname(unlist(dep_list[[sig_celltype]])) & !sim_control %in% unname(unlist(dep_list[[sig_celltype]])))  # sim_celltype and sim_control cannot be dependent on signature cell type

    return(sigs_scores_pvals_tidy_clean)

  }


  sigs_scores_pvals_tidy_clean <- scoreSim(sim, signatures_collection, dep_list, n_null)

  # Get top 25% or top 10 signatures by p-value delta
  top_pval_sigs <- sigs_scores_pvals_tidy_clean %>%
    mutate(pval = -log(pval)) %>%
    mutate(pval_type = ifelse(sig_celltype == sim_celltype, "ctoi", "controls")) %>%
    group_by(sig_celltype, signature, pval_type) %>%
    summarise(pval = mean(pval)) %>%
    pivot_wider(names_from = pval_type, values_from = pval) %>%
    mutate(delta_pval = `ctoi` - `controls`) %>%
    group_by(sig_celltype) %>%
    top_n(n = max(10, 0.25*n()), wt = delta_pval) %>%
    pull(signature)


  # Filter by top_pval_sigs and get top 10% or top 5 signatures by p-value delta
  filtered_sigs <- sigs_scores_pvals_tidy_clean %>%
    filter(signature %in% top_pval_sigs) %>%
    mutate(score_type = ifelse(sig_celltype == sim_celltype, "ctoi", "controls")) %>%
    group_by(sig_celltype, signature, score_type) %>%
    summarise(score = mean(score)) %>%
    pivot_wider(names_from = score_type, values_from = score) %>%
    mutate(delta_score = `ctoi` - `controls`) %>%
    group_by(sig_celltype) %>%
    top_n(n = max(5, 0.1*n()), wt = delta_score) %>%
    pull(signature)


  # celltypes <- unique(gsub("#.*", "", names(signatures_collection)))
  # table(gsub("#.*", "", filtered_sigs))
  # table(gsub("#.*", "", unique(sigs_scores_pvals_tidy_clean$signature)))

  signatures_collection_filtered <- signatures_collection[names(signatures_collection) %in% filtered_sigs]

  out <- list(sim_results = sigs_scores_pvals_tidy_clean,
              filtered_sigs = signatures_collection_filtered)

  return(out)
}
scoreSimulations <- function(signatures, simulations, dep_list, simple){


  celltypes <- names(simulations)
  sims_scored <- pbapply::pblapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]

    if (simple) {
      sim_ranked <- singscore::rankGenes(simulations[[ctoi]])

      scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      }))
      colnames(scores) <- colnames(sim_ranked)

      as_tibble(scores, rownames = "signature") %>%
        pivot_longer(cols = -signature, values_to = "score") %>%
        separate(name, into = c("sim_ct", "sim_frac"), sep = "%%") %>%
        separate(signature, into = "celltype", sep = "#", remove = FALSE, extra = "drop") %>%
        mutate(sim_type = ifelse(celltype == sim_ct, "ctoi", "control"))

    }else{
      lapply(simulations[[ctoi]], function(sim){

        sim_ranked <- singscore::rankGenes(sim)

        scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
          singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
        }))
        colnames(scores) <- colnames(sim)

        as_tibble(scores, rownames = "signature") %>%
          pivot_longer(cols = -signature, values_to = "score") %>%
          separate(name, into = c("sim_ct", "sim_frac"), sep = "%%") %>%
          separate(signature, into = "celltype", sep = "#", remove = FALSE, extra = "drop") %>%
          mutate(sim_type = ifelse(celltype == sim_ct, "ctoi", "control"))

      }) %>%
        bind_rows(., .id = "sim_id")
    }

  })
  names(sims_scored) <- celltypes

  return(sims_scored)

}
getTranformationModels <- function(simulations_scored, RFgamma, XGBparams, modelType){

  set.seed(123)

  fitModel <- function(data, gamma = RFgamma, xgbparams = XGBparams, model_type){


    train_mat <- data %>%
      filter(sim_type == "ctoi") %>%
      select(signature, sim_frac, score) %>%
      mutate(sim_frac = as.numeric(sim_frac)) %>%
      pivot_wider(names_from = signature, values_from = score) %>%
      as.matrix()

    if (model_type == "rf") {
      RF <- RRF::RRF(train_mat[,-1], train_mat[,1], flagReg = 0, importance = TRUE)
      RF_imp <- RF$importance[,"%IncMSE"] / max(RF$importance[,"%IncMSE"])
      RRF <- RRF::RRF(train_mat[,-1], train_mat[,1], flagReg = 1, coefReg = (1-gamma) + gamma*RF_imp)
      return(RRF)
    }


    if (model_type == "xgb") {

      train_mat <- xgboost::xgb.DMatrix(data = train_mat[,-1], label = train_mat[,1])

      # params <- list(
      #   booster = "gbtree",
      #   eta = 0.1,
      #   max_depth = 6,
      #   alpha = 0,  # L1 regularization term
      #   lambda = 1, # L2 regularization term
      #   objective = "reg:logistic"
      # )

      # # Tune parameters
      # cv_results <- xgboost::xgb.cv(params = params, data = train_mat, nrounds = 1000, nfold = 5,
      #                      early_stopping_rounds = 10, verbose = 0)
      # best_nrounds <- cv_results$best_iteration
      #
      # best_params <- list()
      # best_score <- Inf
      # best_reg <- list()
      #
      # for (eta in c(0.01, 0.05)) {
      #   for (max_depth in c(6, 8, 10)) {
      #     for (reg in list(c(0, 0), c(0, 1), c(1, 0), c(1, 1))) {
      #       params$eta <- eta
      #       params$max_depth <- max_depth
      #       params$alpha <- reg[1]
      #       params$lambda <- reg[2]
      #       cv_results <- xgboost::xgb.cv(params = params, data = train_mat, nrounds = best_nrounds, nfold = 5)
      #       mean_test_error <- mean(cv_results$evaluation_log$test_rmse_mean)
      #       if (mean_test_error < best_score) {
      #         best_score <- mean_test_error
      #         best_params <- params
      #         best_reg <- reg
      #       }
      #     }
      #   }
      # }

      model <- xgboost::xgb.train(
        params = xgbparams,
        data = train_mat,
        nrounds = 100
      )

      return(model)
    }

  }

  enframe(simulations_scored, name = "celltype", value = "data") %>%
    rowwise() %>%
    mutate(model = list(fitModel(data, model_type = modelType))) %>%
    dplyr::select(celltype, model) %>%
    return(.)


}
getSpillOverMat <- function(simple_sim, signatures, dep_list, trans_models){

  scoreTransform <- function(mat, signatures, trans_models, is_controls){

    # Score
    mat_ranked <- singscore::rankGenes(mat)
    scores <- t(sapply(signatures, simplify = TRUE, function(sig){
      singscore::simpleScore(mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    }))
    colnames(scores) <- colnames(mat)
    rownames(scores) <- names(signatures)

    transfomed_scores <- as_tibble(scores, rownames = "signatures") %>%
      pivot_longer(cols = -signatures, names_to = "sim_celltype", values_to = "score") %>%
      separate(signatures, into = "sig_celltype", sep = "#", extra = "drop", remove = FALSE) %>%
      group_by(sig_celltype, sim_celltype) %>%
      summarise(mean_score = mean(score)) %>%
      # Transform
      left_join(trans_models, by = c("sig_celltype" = "celltype")) %>%
      rowwise() %>%
      mutate(shifted_score = mean_score - shift_value) %>%
      mutate(transformed_score = round(predict(model, newdata = data.frame("shifted_score" = shifted_score), type = "response"), 2)) %>%
      select(sig_celltype, sim_celltype, transformed_score) %>%
      pivot_wider(names_from = sim_celltype, values_from = transformed_score)

    # mutate(transformed_score = round(predict(model, newdata = data.frame("shifted_score" = shifted_score)) * scaling_value, 2)) %>%
    # select(sig_celltype, mix_ctoi, transformed_score) %>%
    # pivot_wider(names_from = mix_ctoi, values_from = transformed_score)

    # TODO: Inset

    # Convert to matrix
    transfomed_scores_mat <- as.matrix(transfomed_scores[,-1])
    rownames(transfomed_scores_mat) <- pull(transfomed_scores[,1])
    if (!is_controls) {
      transfomed_scores_mat <- transfomed_scores_mat[colnames(mat), colnames(mat)]
    }

    return(transfomed_scores_mat)

  }

  ctoi_mat <- sapply(simple_sim, function(sim){
    sim[,ncol(sim)-1]
  })

  controls_mat <- sapply(simple_sim, function(sim){
    sim[,ncol(sim)]
  })

  colnames(controls_mat) <- unname(sapply(simple_sim, function(sim){
    gsub("%%*.", "", colnames(sim)[ncol(sim)])
  }))

  # Score and transform simulations
  sim_mat_transformed <- scoreTransform(mat = ctoi_mat, signatures, trans_models, is_controls = FALSE)
  controls_mat_uniq <- controls_mat[,!duplicated(colnames(controls_mat))]
  controls_mat_transformed <- scoreTransform(mat = controls_mat_uniq, signatures, trans_models, is_controls = TRUE)
  # Undo unique
  controls_mat_transformed <- sapply(colnames(controls_mat), function(ctrl){
    controls_mat_transformed[,ctrl]
  })
  controls_mat_transformed <- controls_mat_transformed[colnames(sim_mat_transformed), ]

  # Remove control signal from the transformed mixture
  spill_mat <- sim_mat_transformed - controls_mat_transformed

  # Clean and normalize spill matrix
  spill_mat[spill_mat < 0] <- 0
  spill_mat <- spill_mat / diag(spill_mat)

  # Insert zero to dependent cell types
  for(ctoi in rownames(spill_mat)){
    dep_cts <- unname(unlist(dep_list[[ctoi]]))
    dep_cts <- dep_cts[dep_cts != ctoi]
    spill_mat[ctoi, dep_cts] <- 0
  }

  # TODO: Check this parameter
  spill_mat[spill_mat > 0.5] <- 0.5
  diag(spill_mat) <- 1

  # pheatmap::pheatmap(spill_mat, cluster_rows = F, cluster_cols = F)

  return(spill_mat)

}

#' @slot signatures list of xCell2 signatures
#' @slot dependencies list of cell type dependencies
#' @slot transformation_models data frame of cell type transformation models
#' @slot spill_mat matrix of cell types spillover
#' @slot genes_used character vector of genes names used to train the signatures
#' @importFrom methods new
# Create S4 object for the new reference
setClass("xCell2Signatures", slots = list(
  signatures = "list",
  dependencies = "list",
  transformation_models = "data.frame",
  spill_mat = "matrix",
  genes_used = "character"
))


#' xCell2Train function
#'
#' This function generates signatures for each cell type.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import readr
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures VariableFeatures
#' @importFrom Rfast rowMedians rowmeans
#' @importFrom pbapply pblapply pbsapply
#' @importFrom sparseMatrixStats rowMedians
#' @importFrom RRF RRF
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom Matrix rowMeans
#' @importFrom singscore rankGenes simpleScore
#' @param ref A reference gene expression matrix.
#' @param labels A data frame in which the rows correspond to samples in the ref. The data frame must have four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545" or NA if there is no ontology).
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample should match the column name in ref.
#'   "dataset": the cell type sample dataset or subject (for single-cell) as a character.
#' @param data_type Gene expression data type: "rnaseq", "array", or "sc".
#' @param lineage_file Path to the cell type lineage file generated with `xCell2GetLineage` function (optional) .
#' @param sim_fracs A vector of mixture fractions to be used in signature filtering (optional).
#' @param probs A vector of probability thresholds to be used for generating signatures (optional).
#' @param diff_vals A vector of delta values to be used for generating signatures (optional).
#' @param min_genes The minimum number of genes to include in the signature (optional).
#' @param max_genes The maximum number of genes to include in the signature (optional).
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, data_type, lineage_file = NULL, clean_genes = TRUE,
                        sim_fracs = c(0, 0.001, 0.002, 0.004, 0.006, 0.008, seq(0.01, 1, 0.01)), diff_vals = c(1, 1.32, 1.585, 2, 3, 4, 5),
                        min_genes = 5, max_genes = 200, filter_sigs = TRUE, simpleSim = TRUE, sigsFile = NULL, RFgamma = 0.8, XGBparams = list(), modelType = "rf"){


  # Validate inputs
  inputs_validated <- validateInputs(ref, labels, data_type)
  ref <- inputs_validated$ref
  labels <- inputs_validated$labels

  # Prepare reference
  ref <- prepareRef(ref, data_type)

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  pure_ct_mat <- makePureCTMat(ref, labels, use_median = TRUE)
  cor_mat <- getCellTypeCorrelation(pure_ct_mat, data_type)

  # Get cell type dependencies list
  message("Loading dependencies...")
  if (is.null(lineage_file)) {
    dep_list <- xCell2::xCell2GetLineage(labels, out_file = NULL)
  }else{
    dep_list <- getDependencies(lineage_file)
  }

  # Generate/Load signatures
  if (is.null(sigsFile)) {
    # Generate signatures
    if (data_type != "sc") {
      probs <- c(0.01, 0.05, 0.1, 0.25, 0.333, 0.49)
    }else{
      probs <- c(0.1, 0.15, 0.2, 0.25, 0.333, 0.49)
      # Adjust diff values to log1p
      fold_change_vals <- round(2^diff_vals, 4)
      diff_vals <- round(log1p(fold_change_vals - 1), 3)
    }
    message("Calculating quantiles...")
    quantiles_matrix <- makeQuantiles(ref, labels, probs, dep_list, include_descendants = FALSE)
    message("Generating signatures...")
    signatures_collection <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes = TRUE)

    if (!filter_sigs) {
      return(signatures_collection)
    }
  }else{
    message("Loading signatures...")
    signatures_collection <- readRDS(sigsFile)
  }


  # Make simulations
  message("Generating simulations...")
  simple_simulations <- makeSimulations(ref, labels, pure_ct_mat, cor_mat, dep_list, sim_fracs, n_sims = 1, n_samples_sim = NULL, add_noise = FALSE, simple = TRUE)
  if (!simpleSim) {
    complex_simulations <- makeSimulations(ref, labels, pure_ct_mat, cor_mat, dep_list, sim_fracs, n_sims = 5, n_samples_sim = 9, add_noise = TRUE, simple = FALSE)
  }

  # Score simulations
  message("Scoring simulations...")
  if (simpleSim) {
    simple_simulations_scored <- scoreSimulations(signatures = signatures_collection, simulations = simple_simulations, dep_list, simple = TRUE)
  }else{
    complex_simulations_scored <- scoreSimulations(signatures = signatures_collection, simulations = complex_simulations, dep_list, simple = FALSE)
  }

  # TODO: Filter signatures
  message("Filtering signatures...")
  signatures <- signatures_collection


  # Get transformation models
  # TODO: Use filtered signatures
  trans_models <- getTranformationModels(simulations_scored = simple_simulations_scored, RFgamma, XGBparams, modelType)


  # Get spillover matrix
  # TODO: Use filtered signatures
  spill_mat <- matrix()

  # Save results in S4 object
  xCell2Sigs.S4 <- new("xCell2Signatures",
                       signatures = signatures,
                       dependencies = dep_list,
                       transformation_models = trans_models,
                       spill_mat = spill_mat,
                       genes_used = rownames(ref))

  return(xCell2Sigs.S4)

}
