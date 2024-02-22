validateInputs <- function(ref, labels, ref_type){

  if (length(unique(labels$label)) < 3) {
    # TODO: What happens with 2 cell types
    stop("Reference must have at least 3 cell types")
  }

  if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
    stop("ref must be one of those classes: matrix, dgCMatrix, Matrix")
  }

  if (!"data.frame" %in% class(labels)) {
    stop("labels must be a dataframe.")
  }

  if (!ref_type %in% c("rnaseq", "array", "sc")) {
    stop("ref_type should be 'rnaseq', 'array' or 'sc'.")
  }

  if (sum(grepl("_", labels$label)) != 0) {
    message("Changing underscores to dashes in cell-types labels!")
    labels$label <- gsub("_", "-", labels$label)
  }

  out <- list(ref = ref,
              labels = labels)
  return(out)

}
sc2pseudoBulk <- function(ref, labels, min_n_cells, min_ps_samples, seed2use){

  set.seed(seed2use)

  celltypes <- unique(labels$label)

  groups_list <- lapply(celltypes, function(ctoi){

    ctoi_samples <- labels[labels$label == ctoi,]$sample

    # Calculate maximum possible number of groups given min_n_cells
    num_groups <- ceiling(length(ctoi_samples) / min_n_cells)
    if (num_groups < min_ps_samples) {
      num_groups <- min_ps_samples
    }

    # Generate min_ps_samples pseudo samples of CTOI
    if (length(ctoi_samples) > min_ps_samples) {

      ctoi_samples_shuffled <- sample(ctoi_samples, length(ctoi_samples))
      list_of_ctoi_samples_shuffled <- split(ctoi_samples_shuffled, ceiling(seq_along(ctoi_samples_shuffled) / (length(ctoi_samples_shuffled) / num_groups)))

      sapply(list_of_ctoi_samples_shuffled, function(ctoi_group){
        if (length(ctoi_group) == 1) {
          ref[,ctoi_group]
        }else{
          if("matrix" %in% class(ref)) Rfast::rowsums(ref[,ctoi_group]) else Matrix::rowSums(ref[,ctoi_group])
        }
      })

    }else{
      ref[,ctoi_samples]
    }

  })

  names(groups_list) <- celltypes
  pseudo_ref <- as.matrix(bind_cols(groups_list))
  rownames(pseudo_ref) <- rownames(ref)

  pseudo_label <- tibble(labels) %>%
    dplyr::select(ont, label) %>%
    unique() %>%
    right_join(., tibble(label = sub("\\.\\d+$", "", colnames(pseudo_ref)), sample = colnames(pseudo_ref), dataset = "pseudoBulk"), by = "label") %>%
    as.data.frame()

  return(list(ref = pseudo_ref, labels = pseudo_label))

}
normRefMix <- function(ref, mix, filtering_data, ref_type){


  if (ref_type == "sc") {
    message("Normalizing pseudo bulk reference to CPM.")
    # TODO: For non 10X also do TPM
    lib_sizes <- Matrix::colSums(ref)
    norm_factor <- 1e6 / lib_sizes
    ref_norm <- ref %*% Matrix::Diagonal(x = norm_factor)
    colnames(ref_norm) <- colnames(ref)
    ref <- as.matrix(ref_norm) # TODO: Find a way to reduce memory usage by keeping matrix sparse
    message("Scaling reference to log2-space.")
    ref <- log2(ref+1)
  }


  # Select shared genes
  if (is.null(filtering_data)) {
    shared_genes <- intersect(rownames(ref), rownames(mix))
  }else{
    shared_genes <- intersect(intersect(rownames(ref), rownames(mix)), rownames(filtering_data $mixture))
  }
  message(paste0(length(shared_genes), " genes are shared between reference and mixture."))
  ref <- ref[shared_genes,]
  mix <- mix[shared_genes,]


  # Data log-transformation
  if(max(ref) >= 50){
    message("Scaling reference to log2-space.")
    ref <- log2(ref+1)
  }

  if(max(mix) >= 50){
    message("Scaling mixture to log2-space.")
    mix <- log2(mix+1)
  }


  return(list(ref.out = ref, mix.out = mix))
}
makeGEPMat <- function(ref, labels, use_median){

  celltypes <- unique(labels$label)

  gep_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      if(use_median){
        type_vec <- Rfast::rowMedians(as.matrix(ref[,type_samples]))
      }else{
        type_vec <- Rfast::rowmeans(as.matrix(ref[,type_samples]))
      }
    }
  })
  rownames(gep_mat) <- rownames(ref)

  return(gep_mat)
}
getCellTypeCorrelation <- function(gep_mat, ref_type){

  celltypes <- colnames(gep_mat)

  if (ref_type != "sc") {

    # Use top 10% most variable genes
    genes_var <- apply(gep_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.9, na.rm=TRUE)
    gep_mat <- gep_mat[genes_var > most_var_genes_cutoff,]

  }else{

    # Use top 1% most variable genes
    genes_var <- apply(gep_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.99, na.rm=TRUE)
    gep_mat <- gep_mat[genes_var > most_var_genes_cutoff,]

  }


  # Make correlation matrix
  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

  # TODO: Change for loop to apply function to measure time
  for (i in 1:nrow(lower_tri_coord)) {
    celltype_i <- rownames(cor_mat)[lower_tri_coord[i, 1]]
    celltype_j <- colnames(cor_mat)[lower_tri_coord[i, 2]]
    cor_mat[lower_tri_coord[i, 1], lower_tri_coord[i, 2]] <- cor(gep_mat[,celltype_i], gep_mat[,celltype_j], method = "spearman")
    cor_mat[lower_tri_coord[i, 2], lower_tri_coord[i, 1]] <- cor(gep_mat[,celltype_i], gep_mat[,celltype_j], method = "spearman")
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
makeQuantiles <- function(ref, labels, probs, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)
  celltypes <- unique(labels[,2])

  quantiles_mat_list <-  BiocParallel::bplapply(celltypes, function(type){

    type_samples <- labels[,2] == type

    if (sum(type_samples) == 1) {
      # If there is one sample for this cell type -> duplicate the sample to make a data frame
      type.df <- cbind(ref[,type_samples], ref[,type_samples])
    }else{
      type.df <- ref[,type_samples]
    }

    # Calculate quantiles
    # TODO: Balance quantiles by dataset
    type_quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1-probs))), na.rm=TRUE))
  }, BPPARAM = param)
  names(quantiles_mat_list) <- celltypes

  return(quantiles_mat_list)
}
createSignatures <- function(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, top_genes_frac, ncores){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, top_genes_frac){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Set signature thresholds grid
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)
    param.df <- param.df[param.df$diff_vals != 0 | param.df$probs <= 0.05,]
    ntop <- round(length(not_dep_celltypes)*top_genes_frac) # Number of top values given top top_genes_frac of the not dependent cell types
    # Adjusted ntop given the diff value
    param.ranks <- rank((param.df$diff_vals^2)*(1/param.df$probs))
    ntop_adj <- round(1 + ((param.ranks - min(param.ranks)) / (max(param.ranks) - min(param.ranks))) * (ntop - 1)) # Scale param.ranks given ntop
    param.df$ntop <- ntop_adj

    # Generate signatures
    type_sigs <- list()
    for (i in 1:nrow(param.df)){

      # Get a Boolean matrices with genes that pass the quantiles criteria
      diff <- param.df[i, ]$diff_vals # difference threshold
      lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile cutoff


      # Sort upper prob gene value for each not_dep_celltypes
      upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile cutoff
      upper_prob.mat <- sapply(not_dep_celltypes, function(x){
        get(x, quantiles_matrix)[upper_prob,]
      })
      upper_prob.mat <- t(apply(upper_prob.mat, 1, function(x){
        Rfast::Sort(x, descending = TRUE)
      }))

      #  Check diff-prob criteria
      diff_genes.mat <- apply(upper_prob.mat, 2, function(x){
        get(type, quantiles_matrix)[lower_prob,] > x + diff
      })

      # Make signatures
      for (j in 1:param.df[i, ]$ntop) {

        sig_genes <- names(which(diff_genes.mat[,j]))
        n_genes <- length(sig_genes)

        if (n_genes < min_genes) {
          next
        }

        if (n_genes > max_genes) {
          break
        }

        sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, sep = "_")
        type_sigs[[sig_name]] <- sig_genes
      }

    }

    # Remove duplicate signatures
    type_sigs_sorted <- lapply(type_sigs, function(x) sort(x))
    type_sigs_sorted_collapsed <- sapply(type_sigs_sorted, paste, collapse = ",")
    duplicated_sigs <- duplicated(type_sigs_sorted_collapsed)
    type_sigs <- type_sigs[!duplicated_sigs]

    nRelax <- 1
    while (length(type_sigs) < 3 & nRelax <= 5) {
      warning(paste0("Not enough signatures found for ", type, " (relaxing parameters)..."))

      # Relax diff_vals
      param.df$diff_vals <- param.df$diff_vals*0.5

      # Generate signatures
      type_sigs <- list()
      for (i in 1:nrow(param.df)){

        # Get a Boolean matrices with genes that pass the quantiles criteria
        diff <- param.df[i, ]$diff_vals # difference threshold
        lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile cutoff


        # Sort upper prob gene value for each not_dep_celltypes
        upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile cutoff
        upper_prob.mat <- sapply(not_dep_celltypes, function(x){
          get(x, quantiles_matrix)[upper_prob,]
        })
        upper_prob.mat <- t(apply(upper_prob.mat, 1, function(x){
          Rfast::Sort(x, descending = TRUE)
        }))

        #  Check diff-prob creteria
        diff_genes.mat <- apply(upper_prob.mat, 2, function(x){
          get(type, quantiles_matrix)[lower_prob,] > x + diff
        })

        # Make signatures
        for (j in 1:param.df[i, ]$ntop) {

          sig_genes <- names(which(diff_genes.mat[,j]))
          n_genes <- length(sig_genes)

          if (n_genes < min_genes) {
            next
          }

          if (n_genes > max_genes) {
            break
          }

          sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, sep = "_")
          type_sigs[[sig_name]] <- sig_genes
        }

      }
      type_sigs_sorted <- lapply(type_sigs, function(x) sort(x))
      type_sigs_sorted_collapsed <- sapply(type_sigs_sorted, paste, collapse = ",")
      duplicated_sigs <- duplicated(type_sigs_sorted_collapsed)
      type_sigs <- type_sigs[!duplicated_sigs]

      nRelax <- nRelax + 1

    }

    if (length(type_sigs) < 3) {
      errorCondition(paste0("Error: Not enough signatures found for ", type, "!"))
    }

    return(type_sigs)
  }

  param <- BiocParallel::MulticoreParam(workers = ncores)
  celltypes <- unique(labels[,2])

  all_sigs <- BiocParallel::bplapply(celltypes, function(type){

    type.sigs <- getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, top_genes_frac)

    # Check for minimum 3 signatures per cell type
    if (length(type.sigs) < 3) {
      # Relax parameters
      top_genes_frac.tmp <- top_genes_frac
      diff_vals.tmp <- diff_vals
      for (relax_frac in c(0.9, 0.8, 0.7, 0.6, 0.5)) {
        top_genes_frac.tmp <- top_genes_frac.tmp+(1-relax_frac)
        diff_vals.tmp <- diff_vals.tmp*relax_frac
        type.sigs <- getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals.tmp, min_genes, max_genes, top_genes_frac)
        if (length(type.sigs) >= 3) {
          break
        }
      }
    }

    type.sigs
  }, BPPARAM = param)


  all_sigs <- unlist(all_sigs, recursive = FALSE)


  if (length(all_sigs) == 0) {
    stop("No signatures found for reference!")
  }


  return(all_sigs)
}
addEssentialGenes <- function(ref, signatures){

  data("celltype.data", package = "xCell2")
  celltypes <- unique(gsub("#.*", "", names(signatures)))

  for (ct in celltypes) {

    # Check if cell type exists in celltype.data
    if (!ct %in% celltype.data$all_labels) {
      next
    }

    # Get essential genes
    ct_label <- celltype.data[celltype.data$all_labels == ct,]$xCell2_labels
    essen_genes <- unique(unlist(celltype.data[celltype.data$xCell2_labels == ct_label,]$essential_genes))

    if (all(is.na(essen_genes))) {
      next
    }

    essen_genes <- intersect(rownames(ref), essen_genes) # Use only essential genes which are in ref

    # Add to signature
    ct_sigs <- which(startsWith(names(signatures), paste0(ct, "#")))
    for (sig in ct_sigs) {
      signatures[sig][[1]] <- unique(c(signatures[sig][[1]], essen_genes))
    }

  }

  return(signatures)

}
weightFiltData <- function(mix, filtering_data, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)

  celltypes <- unique(labels$label)
  shared_cts <- intersect(celltypes, rownames(filtering_data$truth))
  filtering_mat <- filtering_data$mixture

  # Filtering data log-transformation
  if(max(filtering_mat) >= 50){
    filtering_mat <- log2(filtering_mat+1)
  }

  ds_simi_list <- BiocParallel::bplapply(shared_cts, function(ctoi){

    ctoi_samples <- colnames(filtering_data$truth)[!is.na(filtering_data$truth[ctoi,])]
    ctoi_filt_data <- filtering_mat[,ctoi_samples]
    filt_ds <- unique(gsub("#.*", "", colnames(ctoi_filt_data)))
    if (length(filt_ds) == 1) {
      return(tibble(celltype = ctoi, ds = filt_ds, similarity_scaled = 1))
    }


    # Combine filtering datasets with mixtures
    genes2use <- intersect(rownames(ref), rownames(ctoi_filt_data))
    mix_tmp <- mix
    colnames(mix_tmp) <- paste0("mix#", 1:ncol(mix_tmp))
    filtMix  <- cbind(ctoi_filt_data[genes2use,], mix_tmp[genes2use,])
    datasets <- gsub("#.*", "", colnames(filtMix))
    if(max(mix) >= 50){
      topVarGenes <- names(sort(apply(mix, 1, var), decreasing = TRUE)[1:5000])
    }else{
      topVarGenes <- names(sort(apply(2^mix-1, 1, var), decreasing = TRUE)[1:5000])
    }
    filtMix <- filtMix[topVarGenes,]

    # Run PCA
    pcaResults <- prcomp(t(filtMix))
    pc1Var <- round(summary(pcaResults)$importance[1], 2)
    pc2Var <- round(summary(pcaResults)$importance[2], 2)
    # qplot(pcaResults$x[,1], pcaResults$x[,2], col=datasets, size=4) +
    # labs(x=paste0("PC-1 (", pc1Var, "%)"),
    #      y=paste0("PC-2 (", pc2Var, "%)"))
    pca_scores <- pcaResults$x[,1:2]

    ds_pc_coord <- as_tibble(pca_scores) %>%
      mutate(ds = datasets) %>%
      pivot_longer(-ds, names_to = "PC", values_to = "score") %>%
      group_by(PC, ds) %>%
      summarise(mean_score = mean(score), .groups = "drop_last") %>%
      group_by(PC) %>%
      mutate(mean_score = mean_score + abs(min(mean_score))) %>%
      ungroup()

    ctoi_ds_simi <- ds_pc_coord %>%
      rowwise() %>%
      mutate(dist = abs(ds_pc_coord[ds_pc_coord$ds == "mix" & ds_pc_coord$PC == PC,]$mean_score - mean_score)) %>%
      mutate(dist = ifelse(PC == "PC1", dist*(pc1Var/100),
                           dist*(pc2Var/100))) %>%
      group_by(ds) %>%
      summarise(similarity = (1/sum(dist))*100) %>%
      filter(ds != "mix") %>%
      arrange(-similarity) %>%
      mutate(celltype = ctoi) %>%
      select(celltype, ds, similarity) %>%
      mutate(similarity_scaled = similarity/max(similarity))

    ctoi_ds_simi

  }, BPPARAM = param)
  names(ds_simi_list) <- shared_cts

  filt_ds_weighted <- bind_rows(ds_simi_list)
  return(filt_ds_weighted)

}
filterSignatures <- function(ref, labels, filtering_data, ds_weighted, signatures, top_sigs_frac, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)


  celltypes <- unique(labels$label)
  shared_cts <- intersect(celltypes, rownames(filtering_data$truth))

  filt_sigs <- BiocParallel::bplapply(shared_cts, function(ctoi){

    # Get CTOI signatures
    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

    # Use only samples with truth values
    ctoi_samples <- colnames(filtering_data$truth)[!is.na(filtering_data$truth[ctoi,])]

    # Get filtering data
    ctoi_filt_data <- filtering_data$mixture[, ctoi_samples]

    # Use shared genes
    genes2use <- intersect(rownames(ref), rownames(ctoi_filt_data))
    ctoi_filt_data <- ctoi_filt_data[genes2use,]

    # Calculate CTOI correlations for each dataset in the filtering data
    ctoi_datasets <- gsub("#.*", "", ctoi_samples)
    ds_cors_list <- lapply(unique(ctoi_datasets), function(ds){

      # Get dataset's samples
      ctoi_filt_data_ds <- ctoi_filt_data[,ctoi_datasets == ds]

      #  Rank dataset
      mix_ranked <- singscore::rankGenes(ctoi_filt_data_ds)

      # Score
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){

        suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
      })

      # Calculate correlations
      fracs <- filtering_data$truth[ctoi, colnames(mix_ranked)]
      cors <- apply(scores, 2, function(x){
        cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
      })

      return(sort(cors, decreasing = T))

    })
    names(ds_cors_list) <- unique(ctoi_datasets)

    # Fisher transformation
    z_values <- enframe(ds_cors_list, name = "ds") %>%
      unnest_longer(value, values_to = "rho", indices_to = "sig") %>%
      mutate(z = 0.5 * log((1 + rho) / (1 - rho))) %>%
      mutate(celltype = ctoi)

    # Filter signatures
    sigs_scored <- z_values %>%
      left_join(ds_weighted, by = c("ds", "celltype")) %>%
      mutate(sig_score = z * similarity_scaled) %>%
      group_by(sig) %>%
      summarise(sig_score = mean(sig_score))

    best_sigs <- sigs_scored %>%
      top_frac(top_sigs_frac, wt = sig_score) %>%
      pull(sig)

    if (length(best_sigs) < 10) {
      best_sigs <- sigs_scored %>%
        top_n(10, wt = sig_score) %>%
        pull(sig)
    }


    return(best_sigs)

  }, BPPARAM = param)
  names(filt_sigs) <- shared_cts


  for(ctoi in shared_cts){
    ctoi_sigs <- names(signatures)[startsWith(names(signatures), paste0(ctoi, "#"))]
    sigs2remove <- ctoi_sigs[!ctoi_sigs %in% filt_sigs[[ctoi]]]
    signatures <- signatures[!names(signatures) %in% sigs2remove]
  }

  message("> Signatures from ", length(shared_cts), " cell types have been filtered.")
  out <- list(filt_sigs = signatures,
              filt_cts = shared_cts)

  return(out)
}
makeSimulations <- function(ref, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims, ncores, noise_level = NULL, seed2use){

  set.seed(seed2use)
  param <- BiocParallel::MulticoreParam(workers = ncores)

  celltypes <- unique(labels$label)

  gep_mat_linear <- 2^gep_mat
  if (round(min(gep_mat_linear)) == 1) {
    gep_mat_linear <- gep_mat_linear-1
  }


  sim_list <- BiocParallel::bplapply(celltypes, function(ctoi){

    # Generate CTOI fractions matrix
    ctoi_mat <- matrix(rep(gep_mat_linear[,ctoi], length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs), dimnames = list(rownames(gep_mat_linear), sim_fracs))
    ctoi_mat_frac <- ctoi_mat %*% diag(sim_fracs)

    # Get control cell types
    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
    controls <- celltypes[!celltypes %in% dep_cts]
    if (length(controls) == 0) {
      controls <- names(sort(cor_mat[ctoi,])[1])
    }


    # Generate n_sims simulations
    ctoi_sim_list <- lapply(1:n_sims, function(i){

      # Sample 3-9 non dependent control cell types (if possible)
      if (length(controls) > 3) {
        controls2use <- sample(controls, sample(3:min(8, length(controls)), 1), replace = FALSE)
      }else{
        controls2use <- controls
      }

      # Add one dependent control cell type or closely related cell type to mimic some spillover effect
      if (length(dep_cts[dep_cts != ctoi]) > 0 ) {
        dep_control <- names(sort(cor_mat[ctoi, dep_cts[dep_cts != ctoi]], decreasing = FALSE)[1])
        controls2use <- c(controls2use, dep_control)
      }else{
        dep_control <- names(sort(cor_mat[ctoi,], decreasing = TRUE)[2])
        controls2use <- unique(c(controls2use, dep_control))
      }

      # Generate controls matrix
      if (length(controls2use) > 1) {
        controls_mat <- gep_mat_linear[,controls2use]
        numbers <- runif(ncol(controls_mat)) # Get random numbers for fractions
        controls_mat_frac <- sapply(1-sim_fracs, function(s){
          fracs <- numbers / sum(numbers) * s
          rowSums(controls_mat %*% diag(fracs))
        })
      }else{
        controls_mat <- matrix(rep(gep_mat_linear[,controls2use], length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs), dimnames = list(rownames(gep_mat_linear), sim_fracs))
        controls_mat_frac <- controls_mat %*% diag(1-sim_fracs)
      }


      # Combine fractions
      simulation <- ctoi_mat_frac + controls_mat_frac
      colnames(simulation) <- paste0("mix", "%%", sim_fracs)


      # Add noise
      if (!is.null(noise_level)) {

        # y <- cbind(log2(simulation+1), mix)
        # y <- limma::normalizeBetweenArrays(y)
        # simulation <- y[,1:ncol(simulation)]

        # eps_sigma <- sd(as.vector(simulation)) * noise_level
        eps_sigma <- sd(as.vector(mix)) * noise_level
        eps_gaussian <- array(rnorm(prod(dim(simulation)), 0, eps_sigma), dim(simulation))
        simulation_noised <- 2^(log2(simulation+1) + eps_gaussian)
        colnames(simulation_noised) <- paste0("mix", "%%", sim_fracs)

        #simulation_noised <- 2^(simulation + eps_gaussian)-1

        return(simulation_noised)
      }else{
        simulation
      }

    })
    return(do.call(cbind, ctoi_sim_list))

  }, BPPARAM = param)
  names(sim_list) <- celltypes

  return(sim_list)

}
scoreSimulations <- function(signatures, simulations, n_sims, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)
  celltypes <- names(simulations)
  normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }

  sims_scored <- BiocParallel::bplapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]
    ctoi_sim <- simulations[[ctoi]]
    sim_ranked <- singscore::rankGenes(ctoi_sim)
    colnames(sim_ranked) <- make.unique(colnames(sim_ranked), sep = "%")

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })


    # Normalize scores
    indices <- cut(seq_along(scores[,1]), breaks=n_sims, labels = FALSE)
    scores <- apply(scores, 2, function(sample){
      sample_split <- split(sample, indices)
      norm_scores <- as.numeric(sapply(sample_split, function(x){
        normalize(x)
      }))
      return(norm_scores)
    })

    scores <- cbind(scores, frac = as.numeric(gsub("mix%%", "", colnames(ctoi_sim))))


    return(scores)


  }, BPPARAM = param)

  names(sims_scored) <- celltypes

  return(sims_scored)

}
trainModels <- function(simulations_scored, ncores, seed2use){

  set.seed(seed2use)

  fitModel <- function(data){

    predictors <- data[, -ncol(data)]
    response <- data[, ncol(data)]

    # Remove redundant features
    xgb_params <- list(
      booster = "gbtree",
      alpha = 0,          # Lasso
      lambda = 0,            # Ridge
      eta = 0.01,             # Learning rate
      objective = "reg:squarederror",
      max_depth = 6,
      nthread = 1
    )

    model_tmp <- xgboost::xgboost(
      data = predictors,
      label = response,
      params = xgb_params,
      nrounds = 150,
      verbose = 0
    )


    importance_matrix <- as.data.frame(xgboost::xgb.importance(model = model_tmp))

    if (nrow(importance_matrix) > 50) {
      sigs_filt <- importance_matrix[1:50, 1]
    }else{
      sigs_filt <- importance_matrix[,1]
    }

    # Train model
    xgb_params <- list(
      booster = "gbtree",
      alpha = 1,          # Lasso
      lambda = 1,            # Ridge
      eta = 0.01,             # Learning rate
      objective = "reg:squarederror",
      max_depth = 6,
      nthread = 1
    )

    model_final <- xgboost::xgboost(
      data = predictors[,sigs_filt],
      label = response,
      params = xgb_params,
      nrounds = 150,
      verbose = 0
    )


    return(list(model = model_final, sigs_filtered = sigs_filt))

  }


  param <- BiocParallel::MulticoreParam(workers = ncores)

  #start <- Sys.time()
  models_list <- BiocParallel::bplapply(simulations_scored, function(data){
    fitModel(data)
  }, BPPARAM = param)
  #end <- Sys.time()
  #print(end-start)

  # enframe(models_list, name = "celltype") %>%
  #   unnest(value) %>%
  #   return(.)

  enframe(models_list, name = "celltype") %>%
    unnest_longer(value) %>%
    pivot_wider(names_from = value_id, values_from = value) %>%
    return(.)

}
getSpillOverMat <- function(simulations, signatures, dep_list, models, frac2use){

  scoreTransform <- function(mat, signatures, models){

    # Score
    mat_ranked <- singscore::rankGenes(mat)
    scores <- sapply(signatures, simplify = TRUE, function(sig){
      singscore::simpleScore(mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    rownames(scores) <- colnames(mat)


    transfomed_tbl <- models %>%
      rowwise() %>%
      #mutate(predictions = list(round(predict(model, scores[,startsWith(colnames(scores), paste0(celltype, "#"))], s = model$lambda, type = "response")[,1], 4))) %>%
      mutate(predictions = list(round(randomForestSRC::predict.rfsrc(model, newdata = as.data.frame(scores[,startsWith(colnames(scores), paste0(celltype, "#"))]))$predicted, 4))) %>%
      select(celltype, predictions) %>%
      unnest_longer(predictions, indices_to = "sim_celltype") %>%
      pivot_wider(names_from = sim_celltype, values_from = predictions)


    # Convert to matrix
    transfomed_mat <- as.matrix(transfomed_tbl[,-1])
    rownames(transfomed_mat) <- pull(transfomed_tbl[,1])
    colnames(transfomed_mat) <- rownames(transfomed_mat)


    return(transfomed_mat)

  }


  # Get CTOIs matrix with frac2use fraction
  frac_col <- which(endsWith(colnames(simulations[[1]]), paste0("%%", frac2use)))
  ctoi_mat <- sapply(simulations, function(sim){
    if (length(frac_col) > 1) {
      apply(sim[,frac_col], 1, median)
    }else{
      sim[,frac_col]
    }
  })

  # Get control matrix with CTOI fraction = 0
  frac_col <- which(endsWith(colnames(simulations[[1]]), paste0("%%", 0)))
  controls_mat <- sapply(simulations, function(sim){
    if (length(frac_col) > 1) {
      apply(sim[,frac_col], 1, median)
    }else{
      sim[,frac_col]
    }
  })

  # Score and transform simulations
  sim_transformed <- scoreTransform(mat = ctoi_mat, signatures, models)
  controls_mat_transformed <- scoreTransform(mat = controls_mat, signatures, models)

  # Remove control signal from the transformed mixture
  spill_mat <- sim_transformed - controls_mat_transformed

  # Clean and normalize spill matrix
  spill_mat[spill_mat < 0] <- 0
  spill_mat <- spill_mat / diag(spill_mat)
  spill_mat[is.nan(spill_mat)] <- 0


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
#' @slot models data frame of cell type transformation models
#' @slot spill_mat matrix of cell types spillover
#' @slot genes_used character vector of genes names used to train the signatures
#' @importFrom methods new
# Create S4 object for the new reference
setClass("xCell2Signatures", slots = list(
  signatures = "list",
  dependencies = "list",
  models = "data.frame",
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
#' @import BiocParallel
#' @importFrom limma normalizeBetweenArrays
#' @importFrom outliers grubbs.test
#' @importFrom xgboost xgboost xgb.importance
#' @importFrom Rfast rowMedians rowmeans rowsums Sort
#' @importFrom seqgendiff thin_lib thin_diff
#' @importFrom Matrix rowMeans rowSums colSums
#' @importFrom singscore rankGenes simpleScore
#' @param ref A reference gene expression matrix.
#' @param labels A data frame in which the rows correspond to samples in the ref. The data frame must have four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545" or NA if there is no ontology).
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample should match the column name in ref.
#'   "dataset": the cell type sample dataset or subject (for single-cell) as a character.
#' @param ref_type Gene expression data type: "rnaseq", "array", or "sc".
#' @param mix_type Gene expression data type: "rnaseq", "array", or "sc".
#' @param lineage_file Path to the cell type lineage file generated with `xCell2GetLineage` function (optional).
#' @param top_genes_frac description
#' @param sim_fracs A vector of mixture fractions to be used in signature filtering (optional).
#' @param probs A vector of probability thresholds to be used for generating signatures (optional).
#' @param diff_vals A vector of delta values to be used for generating signatures (optional).
#' @param min_genes The minimum number of genes to include in the signature (optional).
#' @param max_genes The maximum number of genes to include in the signature (optional).
#' @param sigsFile description
#' @param return_sigs description
#' @param return_sigs_filt description
#' @param minPBcells description
#' @param minPBgroups description
#' @param ct_sims description
#' @param samples_frac description
#' @param nCores description
#' @param mix description
#' @param simMethod description
#' @param filtering_data description
#' @param top_sigs_frac description
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, mix = NULL, ref_type, filtering_data = NULL, lineage_file = NULL, top_genes_frac = 1, medianGEP = TRUE, seed = 123, probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4),
                        sim_fracs = c(0, seq(0.01, 0.25, 0.01), seq(0.3, 1, 0.05)), diff_vals = round(c(log2(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3),
                        min_genes = 3, max_genes = 150, return_sigs = FALSE, return_sigs_filt = FALSE, sigsFile = NULL, minPBcells = 30, minPBsamples = 10,
                        ct_sims = 50, samples_frac = 0.1, simMethod = "ref_multi", nCores = 1, top_sigs_frac = 0.05){


  # Validate inputs
  inputs_validated <- validateInputs(ref, labels, ref_type)
  ref <- inputs_validated$ref
  labels <- inputs_validated$labels


  # Generate pseudo bulk from scRNA-Seq reference
  if (ref_type == "sc") {
    message("Converting scRNA-seq reference to pseudo bulk...")
    pb_data <- sc2pseudoBulk(ref, labels, min_n_cells = minPBcells, min_ps_samples = minPBsamples, seed2use = seed)
    ref <- pb_data$ref
    labels <- pb_data$labels
  }

  # Normalize reference/mixture
  # *Also use shared genes
  out <- normRefMix(ref, mix, filtering_data, ref_type)
  ref <- out$ref.out
  mix <- out$mix.out

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  gep_mat <- makeGEPMat(ref, labels, use_median = medianGEP)
  cor_mat <- getCellTypeCorrelation(gep_mat, ref_type)

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
    message("Calculating quantiles...")
    quantiles_matrix <- makeQuantiles(ref, labels, probs, ncores = nCores)
    message("Generating signatures...")
    signatures <- createSignatures(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, top_genes_frac, ncores = nCores)
    signatures <- addEssentialGenes(ref, signatures) # Add essential genes

    if (return_sigs) {
      return(signatures)
    }

  }else{
    # Load signatures
    message("Loading signatures...")
    signatures <- readRDS(sigsFile)
  }

  if (!is.null(filtering_data)) {
    message("Weighting external datasets by similarity to mixture input...")
    ds_weighted <- weightFiltData(mix, filtering_data, ncores = nCores)
    message("Filtering signatures by external datasets...")
    out <- filterSignatures(ref, labels, filtering_data, ds_weighted, signatures, top_sigs_frac, ncores = nCores)
  }

  if (return_sigs_filt) {
    return(list(all_sigs = signatures,
                filt_sigs = out$filt_sigs,
                filt_cts = out$filt_cts,
                genes_used = rownames(ref)))
  }

  signatures <- out$filt_sigs


  # Make simulations
  message("Generating simulations...")
  simulations <- makeSimulations(ref, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims = ct_sims, ncores = nCores, noise_level = 0.025, seed2use = seed)


  message("Scoring simulations...")
  simulations_scored <- scoreSimulations(signatures, simulations, n_sims = ct_sims, nCores)


  # Filter signatures and train RF model
  message("Filtering signatures and training models...")
  models <- trainModels(simulations_scored, ncores = nCores, seed2use = seed)
  signatures <- signatures[unlist(models$sigs_filtered)]
  models <- models[,-3]


  # Get spillover matrix
  # message("Generating spillover matrix...")
  # frac2use <- sim_fracs[which.min(abs(sim_fracs - 0.25))]
  # spill_mat <- getSpillOverMat(simulations, signatures_filt, dep_list, models, frac2use)


  # Save results in S4 object
  xCell2Sigs.S4 <- new("xCell2Signatures",
                       signatures = signatures,
                       dependencies = dep_list,
                       models = models,
                       spill_mat = matrix(),
                       genes_used = rownames(ref)) ###################### Important!!!!


  message("Custom xCell2.0 reference ready!")

  return(xCell2Sigs.S4)

}
