validateInputs <- function(ref, labels, data_type){
  if (length(unique(labels$label)) < 3) {
    stop("Reference must have at least 3 different cell types")
  }

  if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
    stop("ref should be as matrix.")
  }

  if (!"data.frame" %in% class(labels)) {
    stop("labels should be as dataframe.")
  }

  if (!data_type %in% c("rnaseq", "array", "sc")) {
    stop("data_type should be rnaseq, array or scrnaseq.")
  }

  if (sum(grepl("_", labels$label)) != 0 | sum(grepl("_", rownames(ref))) != 0) {
    message("Changing underscores to dashes in genes / cell-types labels!")
    labels$label <- gsub("_", "-", labels$label)
    rownames(ref) <- gsub("_", "-", rownames(ref))
  }

  out <- list(ref = ref,
              labels = labels)
  return(out)

}
getTopVariableGenes <- function(ref, min_genes){

  ref.srt <- Seurat::CreateSeuratObject(counts = ref)
  ref.srt <- Seurat::FindVariableFeatures(ref.srt, selection.method = "vst")
  plot1 <- Seurat::VariableFeaturePlot(ref.srt)
  genesVar <- plot1$data$variance.standardized
  names(genesVar) <- rownames(plot1$data)
  genesVar <- sort(genesVar, decreasing = TRUE)

  if (length(genesVar) < min_genes) {
    genesVar <- names(genesVar[1:min_genes])
  }

  return(genesVar)
}
makePureCTMat <- function(ref, labels){

  celltypes <- unique(labels$label)

  pure_ct_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      type_vec <- if("matrix" %in% class(ref)) Rfast::rowMedians(ref[,type_samples]) else sparseMatrixStats::rowMedians(ref[,type_samples])
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

  quantiles_matrix <-  pbapply::pblapply(celltypes, function(type){

    # Include all the descendants of the cell type in the quantiles calculations
    if (include_descendants) {
      descen_cells <- dep_list[[type]]$descendants
      type_samples <- labels[,2] == type | labels[,2] %in% descen_cells
    }else{
      type_samples <- labels[,2] == type
    }

    # If there is one sample for this cell type -> duplicate the sample to make a data frame
    if (sum(type_samples) == 1) {
      type.df <- cbind(ref[,type_samples], ref[,type_samples])
    }else{
      type.df <- ref[,type_samples]
    }

    # Calculate quantiles
    quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1-probs))), na.rm=TRUE))
  })
  names(quantiles_matrix) <- celltypes

  return(quantiles_matrix)
}
createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Signature parameters
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)

    # Generate signatures
    type_sigs <- list()
    for(i in 1:nrow(param.df)){

      # Get a Boolean matrices with genes that pass the quantiles criteria
      diff <- param.df[i, ]$diff_vals # diffrence criteria
      lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile criteria
      upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile criteria

      diff_genes.mat <- sapply(not_dep_celltypes, function(x){
        get(type, quantiles_matrix)[lower_prob,] > get(x, quantiles_matrix)[upper_prob,] + diff
      })


      if(weight_genes){
        # Score genes using weights
        type_weights <- cor_mat[type, not_dep_celltypes]
        gene_scores <- apply(diff_genes.mat, 1, function(x){
          sum(type_weights[which(x)])
        })
      }else{
        gene_scores <- apply(diff_genes.mat, 1, function(x){
          sum(x)
        })
      }


      gene_passed <- gene_scores[gene_scores > 0]

      # If less than min_genes passed move to next parameters
      if (length(gene_passed) < min_genes) {
        next
      }

      # Save signatures
      gene_passed <- sort(gene_passed, decreasing = TRUE)
      gaps <- seq(0, 1, length.out = 40)
      n_sig_genes <- unique(round(min_genes + (gaps ^ 2) * (max_genes - min_genes)))
      for (n_genes in n_sig_genes) {

        if (length(gene_passed) < n_genes) {
          break
        }

        sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, sep = "_")
        type_sigs[[sig_name]] <- GSEABase::GeneSet(names(gene_passed[1:n_genes]), setName = sig_name)
      }

    }

    if (length(type_sigs) == 0) {
      warning(paste0("No signatures found for ", type))
    }

    return(type_sigs)
  }

  celltypes <- unique(labels[,2])

  all_sigs <- pbapply::pblapply(celltypes, function(type){
    getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes)
  })
  all_sigs <- unlist(all_sigs)


  if (length(all_sigs) == 0) {
    warning("No signatures found for reference!")
  }

  # Make GeneSetCollection object
  signatures_collection <- GSEABase::GeneSetCollection(all_sigs)


  return(signatures_collection)
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

  quantiles_matrix <-  pbapply::pblapply(celltypes, function(type){

    # Include all the descendants of the cell type in the quantiles calculations
    if (include_descendants) {
      descen_cells <- dep_list[[type]]$descendants
      type_samples <- labels[,2] == type | labels[,2] %in% descen_cells
    }else{
      type_samples <- labels[,2] == type
    }

    # If there is one sample for this cell type -> duplicate the sample to make a data frame
    if (sum(type_samples) == 1) {
      type.df <- cbind(ref[,type_samples], ref[,type_samples])
    }else{
      type.df <- ref[,type_samples]
    }

    # Calculate quantiles
    quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1-probs))), na.rm=TRUE))
  })
  names(quantiles_matrix) <- celltypes

  return(quantiles_matrix)
}
createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Signature parameters
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)

    # Generate signatures
    type_sigs <- list()
    for(i in 1:nrow(param.df)){

      # Get a Boolean matrices with genes that pass the quantiles criteria
      diff <- param.df[i, ]$diff_vals # diffrence criteria
      lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile criteria
      upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile criteria

      diff_genes.mat <- sapply(not_dep_celltypes, function(x){
        get(type, quantiles_matrix)[lower_prob,] > get(x, quantiles_matrix)[upper_prob,] + diff
      })


      if(weight_genes){
        # Score genes using weights
        type_weights <- cor_mat[type, not_dep_celltypes]
        type_weights[type_weights < 0.001] <- 0.001 # Minimum correlation to fix zero and negative correlations

        gene_scores <- apply(diff_genes.mat, 1, function(x){
          sum(type_weights[which(x)])
        })
      }else{
        gene_scores <- apply(diff_genes.mat, 1, function(x){
          sum(x)
        })
      }


      gene_passed <- gene_scores[gene_scores > 0]

      # If less than min_genes passed move to next parameters
      if (length(gene_passed) < min_genes) {
        next
      }

      # Save signatures
      gene_passed <- sort(gene_passed, decreasing = TRUE)
      gaps <- seq(0, 1, length.out = 40)
      n_sig_genes <- unique(round(min_genes + (gaps ^ 2) * (max_genes - min_genes)))
      for (n_genes in n_sig_genes) {

        if (length(gene_passed) < n_genes) {
          break
        }

        sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, sep = "_")
        type_sigs[[sig_name]] <- GSEABase::GeneSet(names(gene_passed[1:n_genes]), setName = sig_name)
      }

    }

    if (length(type_sigs) == 0) {
      warning(paste0("No signatures found for ", type))
    }

    return(type_sigs)
  }

  celltypes <- unique(labels[,2])

  all_sigs <- pbapply::pblapply(celltypes, function(type){
    getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes)
  })
  all_sigs <- unlist(all_sigs)


  if (length(all_sigs) == 0) {
    warning("No signatures found for reference!")
  }

  # Make GeneSetCollection object
  signatures_collection <- GSEABase::GeneSetCollection(all_sigs)


  return(signatures_collection)
}
filterSignatures <- function(ref, labels, mixture_fractions, dep_list, cor_mat, pure_ct_mat, signatures_collection){

  # Make mixture matrix (from pure_ct_mat) for the first filtering criteria
  makeMixture <- function(pure_ct_mat, cor_mat, dep_list, mix_frac){

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
    controls_abundance <- sort(table(controls), decreasing = TRUE)
    if(controls_abundance[1] == length(celltypes)-1){
      abundant_control <- names(controls_abundance[1])
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

      mixtures_list <- list(mix1 = list(mix_mat = mix_mat, controls_mat = controls_mat),  mix2 = mix_mat2)

    }else{
      mixtures_list <- list(mix1 = list(mix_mat = mix_mat, controls_mat = controls_mat), mix2 = NULL)

    }

    return(mixtures_list)

  }
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


  # (1) First filtering - Top of delta score between CTOI and median score (of all other cell types)

  mix_frac <- mixture_fractions[length(mixture_fractions)] # Use the highest fraction
  mix_list <- makeMixture(pure_ct_mat, cor_mat, dep_list, mix_frac)
  mix_scores_tidy_clean <- scoreMixture(mix_list$mix1$mix_mat, signatures_collection, dep_list)

  median_sig_score <- mix_scores_tidy_clean %>%
    ungroup() %>%
    filter(sig_celltype != mix_celltype) %>%
    group_by(signature) %>%
    summarise(median_score = median(score))

  sigsPassed <- mix_scores_tidy_clean %>%
    ungroup() %>%
    filter(sig_celltype == mix_celltype) %>%
    left_join(median_sig_score, by = "signature") %>%
    drop_na() %>%  # Some CTOI might not have scores with other cell types (because of controls)
    ungroup() %>%
    mutate(delta_score = score - median_score) %>%
    group_by(sig_celltype) %>%
    top_frac(n=0.2, wt=delta_score) %>%
    pull(signature)

  # If a cell types lost in this filtering step because it is the control of all others
  if(!is.null(mix_list$mix2)){

    celltypes <- unique(labels$label)
    lost_ct <- celltypes[!celltypes %in% gsub("#.*", "", sigsPassed)]

    mix_scores_tidy_clean_lost_ct <- scoreMixture(mix_list$mix2, signatures_collection[startsWith(names(signatures_collection), paste0(lost_ct, "#"))], dep_list)

    median_sig_score_lost_ct <- mix_scores_tidy_clean_lost_ct %>%
      ungroup() %>%
      filter(sig_celltype != mix_celltype) %>%
      group_by(signature) %>%
      summarise(median_score = median(score))

    sigsPassed_lost_ct <- mix_scores_tidy_clean_lost_ct %>%
      ungroup() %>%
      filter(sig_celltype == mix_celltype) %>%
      left_join(median_sig_score_lost_ct, by = "signature") %>%
      mutate(delta_score = score - median_score) %>%
      group_by(sig_celltype) %>%
      top_frac(n=0.2, wt=delta_score) %>%
      pull(signature)

    sigsPassed <- c(sigsPassed, sigsPassed_lost_ct)
  }

  # If a cell types lost in this filtering step because every other cell types in dependencies
  celltypes <- unique(labels$label)
  lost_ct <- celltypes[!celltypes %in% gsub("#.*", "", sigsPassed)]
  if(length(lost_ct) > 0){
    lost_sigs <- names(signatures_collection)[startsWith(names(signatures_collection), paste0(lost_ct, "#"))]
    sigsPassed <- c(sigsPassed, lost_sigs)
  }

  # Second filtering - Top Spearman correlation with simulations
  signatures_filtered <- signatures_collection[names(signatures_collection) %in% sigsPassed]

  # Make simulations
  makeSimulations <- function(ref, labels, mixture_fractions, dep_list, cor_mat, n_ct_sim, add_noise, seed = 123){

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
        if(sum(!colnames(cor_mat) %in% dep_cts) == 1){
          control_ct <- colnames(cor_mat)[!colnames(cor_mat) %in% dep_cts]
        }else{
          control_ct <- names(sort(cor_mat[ctoi, !colnames(cor_mat) %in% dep_cts])[1])
        }

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

  sim_list <- makeSimulations(ref, labels, mixture_fractions, dep_list, cor_mat, n_ct_sim = 10, add_noise = FALSE)
  scores_list <- scoreCTOISimulations(signatures = signatures_filtered, sim_list)

  scores_all_sims_tidy <- enframe(scores_list, name = "celltype") %>%
    unnest_longer(value, indices_to = "sim_id", values_to = "scores") %>%
    rowwise() %>%
    mutate(scores = list(mat2tidy(scores))) %>%
    unnest(scores) %>%
    separate(fraction, into = c("control", "fraction"), sep = "%%") %>%
    mutate(fraction = as.numeric(fraction))

  # Filter signatures
  finalSigsPasses <- scores_all_sims_tidy %>%
    group_by(celltype, signature, sim_id) %>%
    summarise(cor = cor(fraction, score, method = "spearman")) %>%
    group_by(celltype, signature) %>%
    summarise(mean_cor = mean(cor)) %>%
    top_frac(n=0.2, wt=mean_cor) %>%
    pull(signature)


  scores_all_sims_tidy_filtered <- scores_all_sims_tidy %>%
    filter(signature %in% finalSigsPasses)


  out <- list(mixtures = mix_list,
              simulations = scores_all_sims_tidy_filtered)

  return(out)
}
getTranformationModels <- function(simulations){


  fitPolyModel <- function(data, celltype, max_degree = 2){

    df <- data.frame(data)

    models <- list()
    for (degree in 1:max_degree) {
      #model <- lm(fraction ~ poly(shifted_score, degree), data = df)
      model <- suppressWarnings(glm(fraction ~ poly(shifted_score, degree), data = df, family = binomial(link = "logit")))
      models[[degree]] <- model
    }

    bic_values <- sapply(models, BIC)
    best_degree <- which.min(bic_values)

    #final_model <- lm(fraction ~ poly(shifted_score, best_degree), data = df)
    final_model <- suppressWarnings(glm(fraction ~ poly(shifted_score, best_degree), data = df, family = binomial(link = "logit")))


    # # Create a sequence of x values for the prediction
    # x_pred <- seq(min(df$shifted_score), max(df$shifted_score), length.out = 100)
    # # Predict y values using the model
    # y_pred <- predict(model, newdata = data.frame(shifted_score = x_pred), type="response")
    # # Plot the original data
    # plot(df$shifted_score, df$fraction, main = paste0("Polynomial Regression ", best_degree ," degree(s) - ", celltype), xlab = "shifted_score", ylab = "fraction", pch = 19)
    # # Add the fitted line
    # lines(x_pred, y_pred, col = "red", lwd = 2)


    return(final_model)
  }


  simulations %>%
    # Median of simulations
    group_by(celltype, signature, fraction) %>%
    summarise(mean_sim_score = median(score)) %>%
    # Mean of signatures
    group_by(celltype, fraction) %>%
    summarise(mean_ct_score = mean(mean_sim_score)) %>%
    # Get shift values
    mutate(shift_value = min(mean_ct_score)) %>%
    # Shift scores
    mutate(shifted_score = mean_ct_score - shift_value) %>%
    select(-mean_ct_score) %>%
    nest(data = c(shifted_score, fraction)) %>%
    rowwise() %>%
    mutate(model = list(fitPolyModel(data, celltype))) %>%
    select(-data) %>%
    return(.)


}
getSpillOverMat <- function(mixtures, signatures_filtered, dep_list, trans_parameters){

  scoreTransform <- function(mat, signatures_filtered, is_controls){

    # Score
    mat_ranked <- singscore::rankGenes(mat)
    scores <- t(sapply(signatures_filtered, simplify = TRUE, function(sig){
      singscore::simpleScore(mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    }))
    colnames(scores) <- colnames(mat)
    rownames(scores) <- names(signatures_filtered)

    transfomed_scores <- as_tibble(scores, rownames = "signatures") %>%
      pivot_longer(cols = -signatures, names_to = "mix_ctoi", values_to = "score") %>%
      separate(signatures, into = "sig_celltype", sep = "#", extra = "drop") %>%
      group_by(sig_celltype, mix_ctoi) %>%
      summarise(mean_score = mean(score)) %>%
      # Transform
      left_join(trans_parameters, by = c("sig_celltype" = "celltype")) %>%
      rowwise() %>%
      mutate(shifted_score = mean_score - shift_value) %>%
      mutate(transformed_score = round(predict(model, newdata = data.frame("shifted_score" = shifted_score), type = "response"), 2)) %>%
      select(sig_celltype, mix_ctoi, transformed_score) %>%
      pivot_wider(names_from = mix_ctoi, values_from = transformed_score)

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

  mix_mat <- mixtures$mix_mat
  colnames(mix_mat) <- gsub("%%.*", "", colnames(mix_mat))
  controls_mat <- mixtures$controls_mat

  # Score and transform mixtures
  mix_mat_transformed <- scoreTransform(mat = mix_mat, signatures_filtered, is_controls = FALSE)
  controls_mat_uniq <- controls_mat[,!duplicated(colnames(controls_mat))]
  controls_mat_transformed <- scoreTransform(mat = controls_mat_uniq, signatures_filtered, is_controls = TRUE)
  # Undo unique
  controls_mat_transformed <- sapply(colnames(controls_mat), function(ctrl){
    controls_mat_transformed[,ctrl]
  })
  controls_mat_transformed <- controls_mat_transformed[colnames(mix_mat), ]

  # Remove control signal from the transformed mixture
  spill_mat <- mix_mat_transformed - controls_mat_transformed

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

#' @slot signatures ...
#' @slot dependencies ...
#' @slot transformation_models ...
#' @slot spill_mat ...
#' @importFrom methods new
# Create S4 object for the new reference
setClass("xCell2Signatures", slots = list(
  signatures = "GeneSetCollection",
  dependencies = "list",
  transformation_models = "data.frame",
  spill_mat = "matrix"
))


#' xCell2Train function
#'
#' This function generates signatures for each cell type.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import Seurat
#' @import Rfast
#' @import pbapply
#' @import sparseMatrixStats
#' @importFrom  GSEABase GeneSetCollection
#' @importFrom  GSEABase GeneSet
#' @import singscore
#' @param ref A reference gene expression matrix.
#' @param labels A data frame in which the rows correspond to samples in the ref. The data frame must have four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545").
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample should match the column name in ref.
#'   "dataset": the cell type sample dataset or subject (for single-cell) as a character.
#' @param data_type Gene expression data type: "rnaseq", "array", or "sc".
#' @param lineage_file (Optional) Path to the cell type lineage file generated with `xCell2GetLineage` function.
#' @param mixture_fractions A vector of mixture fractions to be used in signature filtering.
#' @param probs A vector of probability thresholds to be used for generating signatures.
#' @param diff_vals A vector of delta values to be used for generating signatures.
#' @param min_genes The minimum number of genes to include in the signature.
#' @param max_genes The maximum number of genes to include in the signature.
#' @param return_unfiltered_signatures for development (remove)!
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, data_type, lineage_file = NULL, mixture_fractions = c(0, 0.001, 0.005, seq(0.01, 0.25, 0.03)),
                        probs = seq(0.25, 0.5, 0.05), diff_vals = c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 10, 12, 15, 17),
                        min_genes = 2, max_genes = 200){


  # Validate inputs
  inputs_validated <- validateInputs(ref, labels, data_type)
  ref <- inputs_validated$ref
  labels <- inputs_validated$labels

  # Use only most variable genes for single-cell data
  if (data_type == "sc") {
    message("Looking for most variable genes with Seurat...")
    topVarGenes <- getTopVariableGenes(ref, min_genes = 10000)
    ref <- ref[rownames(ref) %in% topVarGenes,]
  }

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  pure_ct_mat <- makePureCTMat(ref, labels)
  cor_mat <- getCellTypeCorrelation(pure_ct_mat, data_type)

  # Get cell type dependencies list
  message("Loading dependencies...")
  if (is.null(lineage_file)) {
    dep_list <- xCell2::xCell2GetLineage(labels, out_file = NULL)
  }else{
    dep_list <- getDependencies(lineage_file)
  }

  # Generate signatures for each cell type
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs, dep_list, include_descendants = FALSE)
  message("Generating signatures...")
  signatures_collection <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes = TRUE)

  # Filter signatures
  message("Filtering signatures...")
  filterSignatures.out <- filterSignatures(ref, labels, mixture_fractions, dep_list, cor_mat, pure_ct_mat, signatures_collection)
  signatures_filtered <- signatures_collection[names(signatures_collection) %in% filterSignatures.out$simulations$signature]
  mixtures <- filterSignatures.out$mixtures$mix1

  # Get transformation models
  trans_models <- getTranformationModels(simulations = filterSignatures.out$simulations)

  # Get spillover matrix
  spill_mat <- getSpillOverMat(mixtures, signatures_filtered, dep_list, trans_models)

  xCell2Sigs.S4 <- new("xCell2Signatures",
                      signatures = signatures_filtered,
                      dependencies = dep_list,
                      transformation_models = trans_models,
                      spill_mat = spill_mat)

  return(xCell2Sigs.S4)

}
