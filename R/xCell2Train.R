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
normRefMix <- function(ref, mix, ref_type){


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
  shared_genes <- intersect(rownames(ref), rownames(mix))
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
scoreFiltDS <- function(ref, labels, filtering_data, ds_cor_cutoff, signatures, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)

  filt_cts <- intersect(unique(labels$label), unlist(sapply(filtering_data$truth, rownames)))

  if (length(filt_cts) == 0) {
    warningCondition("Not cell types found for filtering signatures.")
    return(NA)
  }

  ct_ds_scores_list <- BiocParallel::bplapply(filt_cts, function(ctoi){
    # Get CTOI signatures
    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

    # Score
    ds_names <- names(filtering_data$mixture)
    ds_scores_list <- lapply(ds_names, function(ds){


      if(!ctoi %in% rownames(filtering_data$truth[[ds]])){
        return(NA)
      }

      ctoi_mix <- filtering_data$mixture[[ds]]
      ctoi_samples <- filtering_data$truth[[ds]][ctoi,]
      ctoi_samples <- names(ctoi_samples[!is.na(ctoi_samples)])
      ctoi_samples <- ctoi_samples[ctoi_samples %in% colnames(ctoi_mix)]
      genes2use <- intersect(rownames(ctoi_mix), rownames(ref))

      #  Rank filtering dataset
      ctoi_filt_ds_ranked <- singscore::rankGenes(ctoi_mix[, ctoi_samples])

      # Score
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        suppressWarnings(singscore::simpleScore(ctoi_filt_ds_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
      })

      # Remove signatures that failed to score
      if (class(scores)[1] == "list") {
        scores <- scores[lengths(scores) != 0]
        scores <- sapply(scores, cbind)
      }

      fracs <- filtering_data$truth[[ds]][ctoi, ctoi_samples]

      return(cbind(scores, "frac" = fracs))

    })
    names(ds_scores_list) <- ds_names
    ds_scores_list <- ds_scores_list[lengths(ds_scores_list) > 1]

    # Get shared signatures
    sigs2use <- Reduce(intersect, lapply(ds_scores_list, colnames))
    ds_scores_list <- lapply(ds_scores_list, function(x){x[,sigs2use]})

    # Calculate max correlations per ds
    ds2maxCor <- lapply(ds_scores_list, function(x){
      max(apply(x[,-ncol(x)], 2, function(y){
        cor(y, x[,ncol(x)], method = "spearman")
      }))
    })

    # Remove ds that max(cor) < ds_cor_cutoff
    ds2use <- names(which(unlist(ds2maxCor) >= ds_cor_cutoff))
    if (length(ds2use) == 0) {
      return(NA)
    }
    ds_scores_list <- ds_scores_list[ds2use]

    return(ds_scores_list)
    #final_training_ds <- Reduce(rbind, ds_scores_list)

  })
  names(ct_ds_scores_list) <- filt_cts


  return(ct_ds_scores_list)

}
filterSignatures <- function(ref, filtDS_scored, signatures, top_sigs_frac, add_essential_genes, ncores){


  param <- BiocParallel::MulticoreParam(workers = ncores)

  filt_cts <- names(filtDS_scored)


  filt_sigs <- BiocParallel::bplapply(filt_cts, function(ctoi){

    if (is.na(filtDS_scored[[ctoi]][1])) {
      return(list(best_sigs = NA,
                  essential_genes = NA))
    }

    # Calculate CTOI correlations for each dataset in the filtering data
    ds_cors_list <- lapply(filtDS_scored[[ctoi]], function(ds){

      # Calculate correlations
      fracs <- ds[,ncol(ds)]
      cors <- apply(ds[,-ncol(ds)], 2, function(x){
        res <- suppressWarnings(cor.test(x, fracs, method = "spearman", use = "pairwise.complete.obs"))
        c(res$estimate, res$p.value)
      })

      return(tibble(sig = colnames(cors), rho = cors[1,], pval = cors[2,]))
    })

    # Get number of samples per dataset
    ds2n_samples <- enframe(sapply(names(filtDS_scored[[ctoi]]), function(ds){
      mix_samples <- colnames(filtering_data$mixture[[ds]])
      truth_samples <- filtering_data$truth[[ds]][ctoi,]
      truth_samples <- names(truth_samples[!is.na(truth_samples)])
      n_samples <- length(intersect(mix_samples, truth_samples))
      n_samples
    }), name = "ds", value = "n_samples")

    ds_sigs_cors <- enframe(ds_cors_list, name = "ds") %>%
      left_join(ds2n_samples, by = "ds") %>%
      unnest(value)


    # Remove non-significant signatures
    ds_sigs_cors <- ds_sigs_cors %>%
      group_by(ds) %>%
      mutate(padj = p.adjust(pval, method = "BH")) %>%
      ungroup() %>%
      filter(padj <= 0.05)

    ds2use <- unique(ds_sigs_cors$ds)


    if (nrow(ds_sigs_cors) == 0) {
      return(list(best_sigs = NA,
                  essential_genes = NA))
    }


    # Find essential genes
    top_sigs <- ds_sigs_cors %>%
      #group_by(ds) %>%
      #top_frac(0.25, wt=rho) %>% # Top 25% correlation per dataset
      group_by(sig) %>%
      summarise(n_sigs = n()) %>%
      mutate(ds_frac = n_sigs/length(ds2use)) %>%
      filter(ds_frac >= 0.5) %>% # Must be in at least 50% of the datasets %>%
      pull(sig) %>%
      unique()


    if (length(top_sigs) == 0) {
      return(list(best_sigs = NA,
                  essential_genes = NA))
    }

    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
    top_genes <- sort(table(unlist(signatures_ctoi[top_sigs])), decreasing = T)/length(top_sigs)
    top_genes <- names(top_genes[top_genes>=0.5]) # Must be in at least 50% of the signatures
    top_genes <- intersect(top_genes, rownames(ref))

    if (length(top_genes) > 0) {

      ds_top_genes_cors <- lapply(ds2use, function(ds){

        true_fracs <- filtering_data$truth[[ds]][ctoi,]
        true_fracs <- true_fracs[!is.na(true_fracs)]
        top_genes_ctoi_mat <- filtering_data$mixture[[ds]]
        ds_top_genes <- intersect(top_genes, rownames(top_genes_ctoi_mat))
        if (length(ds_top_genes) == 0) {
          return(NULL)
        }

        samples <- intersect(names(true_fracs) , colnames(top_genes_ctoi_mat))
        true_fracs <- true_fracs[samples]
        top_genes_ctoi_mat <- top_genes_ctoi_mat[ds_top_genes, samples]

        if (class(top_genes_ctoi_mat)[1] == "numeric") {
          cors <- suppressWarnings(cor.test(top_genes_ctoi_mat, true_fracs, method = "spearman", use = "pairwise.complete.obs"))
          cors <- tibble(genes = names(top_genes_ctoi_mat), rho = cors$estimate, pval = cors$p.value)
        }else{
          cors <- apply(top_genes_ctoi_mat, 1, function(x){
            res <- suppressWarnings(cor.test(x, true_fracs, method = "spearman", use = "pairwise.complete.obs"))
            c(res$estimate, res$p.value)
          })
          cors <- tibble(genes = colnames(cors), rho = cors[1,], pval = cors[2,])
        }

        return(cors)

      })
      names(ds_top_genes_cors) <- ds2use
      ds_top_genes_cors <- ds_top_genes_cors[!sapply(ds_top_genes_cors, function(x){is.null(x[1])})]

      essential_genes <- enframe(ds_top_genes_cors, name = "ds") %>%
        left_join(ds2n_samples, by = "ds") %>%
        unnest(value) %>%
        group_by(ds) %>%
        mutate(padj = p.adjust(pval, method = "BH")) %>%
        ungroup() %>%
        filter(padj <= 0.05) %>%
        mutate(rho = ifelse(rho == 1, 0.99999, rho),
               rho = ifelse(rho == -1, -0.99999, rho)) %>%
        mutate(z = 0.5 * log((1 + rho) / (1 - rho))) %>%
        mutate(weights = log(n_samples)) %>%
        group_by(genes) %>%
        summarise(weighted_z = weighted.mean(x=z, w=weights)) %>%
        mutate(rho_weigted = (exp(2 * weighted_z) - 1) / (exp(2 * weighted_z) + 1)) %>%
        #top_frac(0.5, wt=rho_weigted) %>%
        pull(genes)

      if (length(essential_genes) == 0) {
        essential_genes <- NA
      }

    }else{
      essential_genes <- NA
    }


    # Find best signatures
    rho_weighted_sigs <- ds_sigs_cors[ds_sigs_cors$sig %in% top_sigs,] %>%
      mutate(rho = ifelse(rho == 1, 0.99999, rho),
             rho = ifelse(rho == -1, -0.99999, rho)) %>%
      ungroup() %>%
      mutate(weights = log(n_samples)) %>%
      mutate(z = 0.5 * log((1 + rho) / (1 - rho))) %>%
      group_by(sig) %>%
      summarise(weighted_z = weighted.mean(x=z, w=weights)) %>%
      mutate(rho_weigted = (exp(2 * weighted_z) - 1) / (exp(2 * weighted_z) + 1))


    top_sigs_frac_adjusted <- ifelse(nrow(rho_weighted_sigs)*top_sigs_frac > 10, top_sigs_frac, 10/nrow(rho_weighted_sigs)) # Minimum 10 signatures

    best_sigs <- rho_weighted_sigs %>%
      top_frac(top_sigs_frac_adjusted, wt = rho_weigted) %>%
      pull(sig)

    if(length(best_sigs) == 0){
      best_sigs <- NA
    }



    return(list(best_sigs = best_sigs,
                essential_genes = essential_genes))

  }, BPPARAM = param)
  names(filt_sigs) <- filt_cts


  # Filter genes
  for(ctoi in filt_cts){
    ctoi_best_sigs <- filt_sigs[[ctoi]]$best_sigs
    if (all(!is.na(ctoi_best_sigs))) {
      ctoi_sigs <- names(signatures)[startsWith(names(signatures), paste0(ctoi, "#"))]
      sigs2remove <- ctoi_sigs[!ctoi_sigs %in% ctoi_best_sigs]
      signatures <- signatures[!names(signatures) %in% sigs2remove]
    }

    ctoi_essential_genes <- filt_sigs[[ctoi]]$essential_genes
    if (add_essential_genes & all(!is.na(ctoi_essential_genes))) {
      # Add essential genes
      ctoi_sigs <- names(signatures)[startsWith(names(signatures), paste0(ctoi, "#"))]
      for (sig in ctoi_sigs) {
        signatures[sig][[1]] <- unique(c(signatures[sig][[1]], ctoi_essential_genes))
      }
    }
  }

  # Remove duplicate signatures
  sigs_sorted <- lapply(signatures, function(x) sort(x))
  sigs_sorted_collapsed <- sapply(sigs_sorted, paste, collapse = ",")
  duplicated_sigs <- duplicated(sigs_sorted_collapsed)
  signatures <- signatures[!duplicated_sigs]

  # Number of  cell types passed filtering
  filt_sigs <- filt_sigs[sapply(filt_cts, function(ctoi){
    all(!is.na(filt_sigs[[ctoi]]$best_sigs))
  })]
  filts_ct <- names(filt_sigs)

  message("> Signatures from ", length(filts_ct), " cell types have been filtered.")
  if (add_essential_genes) {
    message("> ", length(unique(unlist(lapply(filt_sigs, function(x){x$essential_genes})))), " essential genes have been added to signatures.")
  }
  out <- list(filt_sigs = signatures,
              filt_cts = filts_ct)

  return(out)
}
filterSignatures2 <- function(ref, labels, filtering_data, signatures, top_sigs_frac, add_essential_genes, ncores){


  param <- BiocParallel::MulticoreParam(workers = ncores)


  celltypes <- unique(labels$label)
  shared_cts <- intersect(celltypes, unlist(sapply(filtering_data$truth, rownames)))


  filt_sigs <- BiocParallel::bplapply(shared_cts, function(ctoi){

    # Get CTOI signatures
    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

    # Calculate CTOI correlations for each dataset in the filtering data
    ds_cors_list <- sapply(names(filtering_data$mixture), function(ds){


      if(!ctoi %in% rownames(filtering_data$truth[[ds]])){
        return(NA)
      }

      ctoi_mix <- filtering_data$mixture[[ds]]
      ctoi_samples <- filtering_data$truth[[ds]][ctoi,]
      ctoi_samples <- names(ctoi_samples[!is.na(ctoi_samples)])
      ctoi_samples <- ctoi_samples[ctoi_samples %in% colnames(ctoi_mix)]
      genes2use <- intersect(rownames(ctoi_mix), rownames(ref))

      #  Rank filtering dataset
      ctoi_filt_ds_ranked <- singscore::rankGenes(ctoi_mix[, ctoi_samples])

      # Score
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        suppressWarnings(singscore::simpleScore(ctoi_filt_ds_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
      })

      # Remove signatures that failed to score
      if (class(scores)[1] == "list") {
        scores <- scores[lengths(scores) != 0]
        scores <- sapply(scores, cbind)
      }

      # Calculate correlations
      fracs <- filtering_data$truth[[ds]][ctoi, ctoi_samples]
      cors <- apply(scores, 2, function(x){
        cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
      })

      return(cors)

    })
    ds_cors_list <- ds_cors_list[lengths(ds_cors_list) > 1] # Remove NAs

    # Get number of samples per dataset
    ds2n_samples <- enframe(sapply(names(ds_cors_list), function(ds){
      mix_samples <- colnames(filtering_data$mixture[[ds]])
      truth_samples <- filtering_data$truth[[ds]][ctoi,]
      truth_samples <- names(truth_samples[!is.na(truth_samples)])
      n_samples <- length(intersect(mix_samples, truth_samples))
      n_samples
    }), name = "ds", value = "n_samples")

    ds_sigs_cors <- enframe(ds_cors_list, name = "ds") %>%
      unnest_longer(value, values_to = "rho", indices_to = "sig")


    # External dataset must max(rho) >= 0.5 to be used in filtering
    ds2use <- ds_sigs_cors %>%
      group_by(ds) %>%
      summarise(max_rho = max(rho)) %>%
      filter(max_rho >= 0.5) %>%
      pull(ds)

    if (length(ds2use) == 0) {
      return(list(best_sigs = NA,
                  essential_genes = NA))
    }

    ds_sigs_cors <- ds_sigs_cors %>%
      filter(ds %in% ds2use)

    # Find essential genes
    top_sigs <- ds_sigs_cors %>%
      group_by(ds) %>%
      top_frac(0.25, wt=rho) %>% # Top 25% correlation per dataset
      #filter(rho >= 0.3) %>%
      group_by(sig) %>%
      summarise(n_sigs = n()) %>%
      mutate(ds_frac = n_sigs/length(ds2use)) %>%
      filter(ds_frac >= 0.5) %>% # Must be in at least 50% of the datasets %>%
      pull(sig) %>%
      unique()

    if (length(top_sigs) == 0) {
      return(list(best_sigs = NA,
                  essential_genes = NA))
    }

    top_genes <- sort(table(unlist(signatures_ctoi[top_sigs])), decreasing = T)/length(top_sigs)
    top_genes <- names(top_genes[top_genes>=0.5]) # Must be in at least 50% of the signatures
    top_genes <- intersect(top_genes, rownames(ref))

    if (length(top_genes) > 0) {
      ds_top_genes_cors <- lapply(ds2use, function(ds){

        true_fracs <- filtering_data$truth[[ds]][ctoi,]
        true_fracs <- true_fracs[!is.na(true_fracs)]
        top_genes_ctoi_mat <- filtering_data$mixture[[ds]]
        ds_top_genes <- intersect(top_genes, rownames(top_genes_ctoi_mat))
        if (length(ds_top_genes) == 0) {
          return(NULL)

        }

        samples <- intersect(names(true_fracs) , colnames(top_genes_ctoi_mat))
        true_fracs <- true_fracs[samples]
        top_genes_ctoi_mat <- top_genes_ctoi_mat[ds_top_genes, samples]

        if (class(top_genes_ctoi_mat)[1] == "numeric") {
          cors <- cor(top_genes_ctoi_mat, true_fracs, method = "spearman", use = "pairwise.complete.obs")
        }else{
          cors <- apply(top_genes_ctoi_mat, 1, function(x){
            cor(x, true_fracs, method = "spearman", use = "pairwise.complete.obs")
          })
        }

        cors[is.na(cors)] <- -1
        return(cors)

      })
      names(ds_top_genes_cors) <- ds2use
      ds_top_genes_cors <- ds_top_genes_cors[!sapply(ds_top_genes_cors, function(x){is.null(x[1])})]

      essential_genes <- enframe(ds_top_genes_cors, name = "ds", value = "rho") %>%
        left_join(ds2n_samples, by = "ds") %>%
        unnest_longer(rho, indices_to = "genes") %>%
        mutate(rho = ifelse(rho == 1, 0.99999, rho),
               rho = ifelse(rho == -1, -0.99999, rho)) %>%
        mutate(z = 0.5 * log((1 + rho) / (1 - rho))) %>%
        mutate(weights = log(n_samples)) %>%
        group_by(genes) %>%
        summarise(weighted_z = weighted.mean(x=z, w=weights)) %>%
        mutate(rho_weigted = (exp(2 * weighted_z) - 1) / (exp(2 * weighted_z) + 1)) %>%
        filter(rho_weigted >= 0.3) %>% # Genes must have at least 0.3 weighted correlation with the filtering data
        top_frac(0.5, wt=rho_weigted) %>%
        pull(genes)

      if (length(essential_genes) == 0) {
        essential_genes <- NA
      }

    }else{
      essential_genes <- NA
    }


    # Find best signatures
    rho_weighted_sigs <- ds_sigs_cors[ds_sigs_cors$sig %in% top_sigs,] %>%
      left_join(ds2n_samples, by = "ds") %>%
      mutate(rho = ifelse(rho == 1, 0.99999, rho),
             rho = ifelse(rho == -1, -0.99999, rho)) %>%
      ungroup() %>%
      mutate(weights = log(n_samples)) %>%
      mutate(z = 0.5 * log((1 + rho) / (1 - rho))) %>%
      group_by(sig) %>%
      summarise(weighted_z = weighted.mean(x=z, w=weights)) %>%
      mutate(rho_weigted = (exp(2 * weighted_z) - 1) / (exp(2 * weighted_z) + 1))

    top_sigs_frac_adjusted <- ifelse(nrow(rho_weighted_sigs)*top_sigs_frac > 10, top_sigs_frac, 10/nrow(rho_weighted_sigs)) # Minimum 10 signatures

    best_sigs <- rho_weighted_sigs %>%
      top_frac(top_sigs_frac_adjusted, wt = rho_weigted) %>%
      pull(sig)

    if(length(best_sigs) == 0){
      best_sigs <- NA
    }

    return(list(best_sigs = best_sigs,
                essential_genes = essential_genes))

  }, BPPARAM = param)
  names(filt_sigs) <- shared_cts


  # Filter genes
  for(ctoi in shared_cts){
    ctoi_best_sigs <- filt_sigs[[ctoi]]$best_sigs
    if (all(!is.na(ctoi_best_sigs))) {
      ctoi_sigs <- names(signatures)[startsWith(names(signatures), paste0(ctoi, "#"))]
      sigs2remove <- ctoi_sigs[!ctoi_sigs %in% ctoi_best_sigs]
      signatures <- signatures[!names(signatures) %in% sigs2remove]
    }

    ctoi_essential_genes <- filt_sigs[[ctoi]]$essential_genes
    if (add_essential_genes & all(!is.na(ctoi_essential_genes))) {
      # Add essential genes
      ctoi_sigs <- names(signatures)[startsWith(names(signatures), paste0(ctoi, "#"))]
      for (sig in ctoi_sigs) {
        signatures[sig][[1]] <- unique(c(signatures[sig][[1]], ctoi_essential_genes))
      }
    }
  }

  # Remove duplicate signatures
  sigs_sorted <- lapply(signatures, function(x) sort(x))
  sigs_sorted_collapsed <- sapply(sigs_sorted, paste, collapse = ",")
  duplicated_sigs <- duplicated(sigs_sorted_collapsed)
  signatures <- signatures[!duplicated_sigs]

  # Number of  cell types passed filtering
  filt_sigs <- filt_sigs[sapply(shared_cts, function(ctoi){
    all(!is.na(filt_sigs[[ctoi]]$best_sigs))
  })]
  filts_ct <- names(filt_sigs)

  message("> Signatures from ", length(filts_ct), " cell types have been filtered.")
  if (add_essential_genes) {
    message("> ", length(unique(unlist(lapply(filt_sigs, function(x){x$essential_genes})))), " essential genes have been added to signatures.")
  }
  out <- list(filt_sigs = signatures,
              filt_cts = filts_ct)

  return(out)
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
makeSimulations <- function(ref, mix, signatures, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims, ncores, noise_level, seed2use){


  sampleControls <- function(ctoi, controls, dep_cts, cor_mat){

    # Sample 3-9 non dependent control cell types (if possible)
    if (length(controls) > 3) {
      controls2use <- sample(controls, sample(3:min(8, length(controls)), 1), replace = FALSE)
    }else{
      controls2use <- controls
    }

    # For half of the simulations add a related cell type to mimic some spillover effect
    if (runif(1) <= 0.5) {
      cts_cors <- cor_mat[ctoi, controls]
      related_cts <- names(which(cts_cors > 0.85 & cts_cors < 0.95))
      if (length(related_cts) > 0) {
        controls2use <- unique(c(controls2use, sample(related_cts, 1)))

      }
    }

    return(controls2use)

  }
  getControls <- function(ctoi, controls, mix_ranked, gep_mat_linear, signatures, dep_cts, sim_fracs, cor_mat, n_sims){


    # Generate sets of different controls combinations
    n_sets <- n_sims
    controls_sets <- lapply(1:n_sets, function(c){
      sampleControls(ctoi, controls, dep_cts, cor_mat)
    })
    names(controls_sets) <- paste0("set-", 1:n_sets)

    # Generate sets of different controls proportions
    controls_props <- lapply(controls_sets, function(p){
      numbers <- runif(length(p))
      fracs2use <- numbers / sum(numbers)
    })


    return(list(c = controls_sets,
                p = controls_props))

  }

  param <- BiocParallel::MulticoreParam(workers = ncores, RNGseed = seed2use)

  celltypes <- unique(labels$label)
  gep_mat_linear <- 2^gep_mat
  if (round(min(gep_mat_linear)) == 1) {
    gep_mat_linear <- gep_mat_linear-1
  }
  mix_ranked <- singscore::rankGenes(mix)

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

    # get set of controls that resemble the mixture's shift value
    controls_sets <- getControls(ctoi, controls, mix_ranked, gep_mat_linear, signatures, dep_cts, sim_fracs, cor_mat, n_sims)

    # Generate n_sims simulations
    ctoi_sim_list <- lapply(1:n_sims, function(i){

      controls2use <- controls_sets$c[[i]]

      if (length(controls2use) > 1) {
        controls_mat <- gep_mat_linear[,controls2use]
      }else{
        controls_mat <- matrix(rep(gep_mat_linear[,controls2use], length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs), dimnames = list(rownames(gep_mat_linear), sim_fracs))
      }

      fracs2use <- controls_sets$p[[i]]

      controls_mat_frac <- sapply(1-sim_fracs, function(s){
        perturbed_vec <- -1
        while (min(perturbed_vec) < 0) {
          perturbed_vec <- fracs2use + runif(length(fracs2use), min=-0.1, max=0.1) # Change the proportions a little bit

        }
        fracs <- s * perturbed_vec / sum(perturbed_vec)
        rowSums(controls_mat %*% diag(fracs))
      })


      # Combine fractions
      simulation <- ctoi_mat_frac + controls_mat_frac
      colnames(simulation) <- paste0("mix", "%%", sim_fracs)


      # Add noise
      if (!is.null(noise_level)) {

        eps_sigma <- sd(as.vector(mix)) * noise_level
        eps_gaussian <- array(rnorm(prod(dim(simulation)), 0, eps_sigma), dim(simulation))
        simulation_noised <- 2^(log2(simulation+1) + eps_gaussian)
        colnames(simulation_noised) <- paste0("mix", "%%", sim_fracs)

        return(simulation_noised)
      }else{
        simulation
      }

    })
    names(ctoi_sim_list) <- paste0("sim-", 1:length(ctoi_sim_list))
    ctoi_sim_list

  }, BPPARAM = param)
  names(sim_list) <- celltypes

  return(sim_list)

}
scoreSimulations_old <- function(signatures, simulations, n_sims, sim_fracs, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)
  celltypes <- names(simulations)

  sims_scored <- BiocParallel::bplapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]
    ctoi_sim_list <- simulations[[ctoi]]

    ctoi_sims <- lapply(1:length(ctoi_sim_list), function(sim_id){

      sim <- ctoi_sim_list[[sim_id]]
      sim_ranked <- singscore::rankGenes(sim)
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      })
      scores <- cbind(scores, frac = sim_fracs)

      as_tibble(scores) %>%
        pivot_longer(cols = -frac, names_to = "sig", values_to = "score") %>%
        mutate(sim = names(ctoi_sim_list)[sim_id])

    })
    ctoi_sims <- bind_rows(ctoi_sims)

    ctoi_sims %>%
      mutate(celltype = ctoi)

  }, BPPARAM = param)
  names(sims_scored) <- celltypes


  return(sims_scored)

}
scoreSimulations <- function(signatures, simulations, n_sims, sim_fracs, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)
  celltypes <- names(simulations)

  sims_scored <- BiocParallel::bplapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]
    ctoi_sim_list <- simulations[[ctoi]]

    ctoi_sims <- lapply(1:length(ctoi_sim_list), function(sim_id){

      sim <- ctoi_sim_list[[sim_id]]
      sim_ranked <- singscore::rankGenes(sim)
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      })
      scores <- cbind(scores, frac = sim_fracs)

      scores

    })
    ctoi_sims <- Reduce(rbind, ctoi_sims)

    ctoi_sims
  }, BPPARAM = param)
  names(sims_scored) <- celltypes


  return(sims_scored)

}
learnParams <- function(simulations_scored, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)

  celltypes <- names(simulations_scored)

  linearParams <- BiocParallel::bplapply(celltypes, function(ctoi){
    tp <- simulations_scored[[ctoi]] %>%
      group_by(celltype, sig, frac) %>%
      summarise(score = mean(score), .groups = "drop") %>%
      group_by(celltype, frac) %>%
      summarise(score = mean(score), .groups = "drop") %>%
      mutate(score = score - min(score)) %>%
      filter(frac > 0) %>%
      summarise(tp = list(try(minpack.lm::nlsLM(score ~ a * frac^b, start = list(a=1, b=1), control = list(maxiter = 500)), silent = TRUE))) %>%
      pull(tp)

    return(tp[[1]])
  }, BPPARAM = param)
  names(linearParams) <- celltypes

  enframe(linearParams, name = "celltype") %>%
    rowwise() %>%
    mutate(a = coef(value)[[1]],
           b = coef(value)[[2]]) %>%
    select(-value) %>%
    ungroup() %>%
    return(.)

}
linearTransform <- function(params, simulations_scored, filtering_data, signatures, min_filt_ds, ref, ncores){


  # Those will be use for fitting linear model and spillover correction
  simulations_transformed <- simulations_scored %>%
    bind_rows() %>%
    group_by(celltype, sig, frac) %>%
    summarise(score = mean(score), .groups = "drop") %>%
    group_by(celltype, frac) %>%
    summarise(score = mean(score), .groups = "drop_last") %>%
    mutate(score = score - min(score)) %>%
    left_join(params, by = "celltype") %>%
    group_by(celltype) %>%
    mutate(score = (score^(1/b)) / a) %>%
    select(-c(a, b)) %>%
    pivot_wider(names_from = celltype, values_from = score) %>%
    as.data.frame()

  # Those will be use for fitting multiple linear regression model for the filtered cell types
  if (!is.null(filtering_data)) {
    param <- BiocParallel::MulticoreParam(workers = ncores)

    filt_cts <- intersect(names(simulations_scored), unlist(sapply(filtering_data$truth, rownames)))
    filtering_transformed_list <- BiocParallel::bplapply(filt_cts, function(ctoi){

      signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

      ds2use <- sapply(filtering_data$truth, function(ds){
        ctoi %in% rownames(ds)
      })

      if (sum(ds2use) < min_filt_ds) {
        return(NULL)
      }

      ds_scores_list <- lapply(names(which(ds2use)), function(ds){

        ctoi_mix <- filtering_data$mixture[[ds]]
        ctoi_samples <- filtering_data$truth[[ds]][ctoi,]
        ctoi_samples <- names(ctoi_samples[!is.na(ctoi_samples)])
        ctoi_samples <- ctoi_samples[ctoi_samples %in% colnames(ctoi_mix)]
        genes2use <- intersect(rownames(ctoi_mix), rownames(ref))

        #  Rank filtering dataset
        ctoi_filt_ds_ranked <- singscore::rankGenes(ctoi_mix[genes2use, ctoi_samples])


        # Score
        scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
          suppressWarnings(singscore::simpleScore(ctoi_filt_ds_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
        })

        # Remove signatures that failed to score
        if (class(scores)[1] == "list") {
          scores <- scores[lengths(scores) != 0]
          scores <- sapply(scores, cbind)
        }

        fracs <- filtering_data$truth[[ds]][ctoi, ctoi_samples]

        cbind(frac = fracs, scores)

      })
      names(ds_scores_list) <- names(which(ds2use))

      # Balance the training data

      # Sort datasets by number of samples
      nsamples <- sapply(ds_scores_list, nrow)
      ds_scores_list <- ds_scores_list[order(nsamples)]
      min_size <- min(nsamples)

      # Split true fracs to bins of min_size samples
      fracs_list <- sapply(ds_scores_list, function(ds){ds[,1]})
      fracs <- unlist(fracs_list)

      q_bins <- cut(unlist(fracs), breaks = quantile(fracs, probs = seq(0, 1, length.out = min_size + 1)), include.lowest = TRUE)
      names(q_bins) <- names(unlist(fracs))

      # Balance datasets
      samples2use <- c()
      for (q in unique(q_bins)) {

        qbinsampsle2use <- c()
        for (ds in fracs_list) {
          qbinsampsle2use <- c(qbinsampsle2use, names(ds)[names(ds) %in% names(q_bins[q_bins == q])])
          if (length(qbinsampsle2use) >=  min_size) {
            qbinsampsle2use <- qbinsampsle2use[1:min_size]
            break
          }
        }
        samples2use <- c(samples2use, qbinsampsle2use)

      }


      sigs2use <- Reduce(intersect, lapply(ds_scores_list, colnames))
      ds_scores_list <- lapply(ds_scores_list, function(x){x[,sigs2use]})
      ds_scores_mat <- Reduce(rbind, ds_scores_list)

      fracs <- ds_scores_mat[samples2use, 1]
      # Range of fracs must be at least 3% and at least 10 samples
      if (max(fracs) - min(fracs) < 0.03 | length(fracs) < 10) {
        return(NULL)
      }
      scores <- ds_scores_mat[samples2use, -1]

      # Transform
      a <- params[params$celltype == ctoi,]$a
      b <- params[params$celltype == ctoi,]$b

      scores_transformed <- apply(scores, 2, function(x){
        x_shifted <- x-min(x)
        (x_shifted^(1/b)) / a
      })

      df <- as.data.frame(cbind(fracs = fracs, scores_transformed))
      return(df)

    }, BPPARAM = param)
    names(filtering_transformed_list) <- filt_cts
    filtering_transformed_list <- Filter(Negate(is.null), filtering_transformed_list)

  }else{
    filtering_transformed_list <- NA
  }

  return(list(sims = simulations_transformed,
              filt = filtering_transformed_list))

}
trainModels_old <- function(data_transfomed, params, ncores){


  fitModel <- function(data){

    if (ncol(data) == 2) {
      # Linear model from simulations
      names(data)[2] <- "score"
      model <- lm(frac~., data)
      return(model)
    }else{
      # Ridge model from filtering data
      model <- glmnet::cv.glmnet(x = as.matrix(data[,-1]), y = data[,1], alpha=0)
      return(model)
    }

  }


  param <- BiocParallel::MulticoreParam(workers = ncores)
  celltypes <- colnames(data_transfomed$sims)[-1]

  models_list <- BiocParallel::bplapply(celltypes, function(ctoi){

    if (ctoi %in%  names(data_transfomed$filt)) {
      data <- data_transfomed$filt[[ctoi]]
    }else{
      data <- data_transfomed$sims[,c("frac", ctoi)]
    }

    fitModel(data)
  }, BPPARAM = param)
  names(models_list) <- celltypes

  enframe(models_list, name = "celltype", value = "model") %>%
    return(.)

}
trainModels <- function(simulations_scored, filtering_data = NA, signatures, L2, ncores){

  fitModel <- function(data, nCores = ncores, l2 = L2){


    # Fix for special JSON characters
    feature_names_filt <- colnames(data)
    special_json <- grepl(",", feature_names_filt)
    if (any(special_json)) {
      colnames(data) <- gsub(pattern = ",", "", colnames(data))
    }

    predictors <- data[, -ncol(data)]
    response <- data[, ncol(data)]
    predictors <- apply(predictors, 2, function(x){x-min(x)})
    predictors <- apply(predictors, 2, scale)


    # Tune
    # Determine dynamic hyperparameter ranges based on data size
    num_samples <- nrow(data)
    max_leaves <- max(15, min(255, round(sqrt(num_samples) * 2)))
    min_data_range <- round(seq(from = max(1, num_samples / 1000), to = max(5, num_samples / 50), length.out = 3))

    # Define the dynamic range of hyperparameters to tune
    tuning_grid <- expand.grid(
      num_leaves = round(seq(from = max(15, max_leaves / 2), to = max_leaves, length.out = 3)),
      min_data_in_leaf = min_data_range,
      learning_rate = c(0.005, 0.01, 0.05, 0.1),
      lambdasL2 = c(0, 0.1, 0.5, 1, 10, 100),
      stringsAsFactors = FALSE
    )

    nfold <- min(max(5, round(num_samples / 100)), 10)

    best_score <- Inf
    best_params <- NULL

    # Tune hyperparameters
    for (i in 1:nrow(tuning_grid)) {
      params <- list(
        objective = "regression",
        metric = "rmse",
        num_threads = nCores,
        learning_rate = tuning_grid$learning_rate[i],
        num_leaves = tuning_grid$num_leaves[i],
        min_data_in_leaf = tuning_grid$min_data_in_leaf[i],
        lambda_l2 = tuning_grid$lambdasL2[i],
        feature_pre_filter = FALSE
      )

      dtrain <- lightgbm::lgb.Dataset(data = predictors, label = response, params = params)


      cv_result <- lightgbm::lgb.cv(
        params = params,
        data = dtrain,
        nfold = nfold,
        nrounds = 500,
        early_stopping_rounds = 10,
        verbose = -1
      )

      # Update best score and parameters if improvement is found
      if (cv_result$best_score < best_score) {
        best_score <- cv_result$best_score
        best_params <- params
      }
    }


    dtrain <- lightgbm::lgb.Dataset(data = predictors, label = response, params = best_params)

    model <- lightgbm::lgb.train(
      params = best_params,
      data = dtrain,
      nrounds = 1000
    )


    return(model)

  }


  celltypes <- names(simulations_scored)


  models_list <- lapply(celltypes, function(ctoi){


    if (ctoi %in% filt_cts) {
      data <- Reduce(rbind, filtDS_scored_filt[[ctoi]][ds2use])
    }else{
      data <- simulations_scored[[ctoi]]
    }

    ctoi_model <- fitModel(data)

    return(tibble(celltype = ctoi, model = list(ctoi_model)))

  })
  names(models_list) <- celltypes

  bind_rows(models_list) %>%
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
  params = "data.frame",
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
#' @importFrom minpack.lm nlsLM
#' @importFrom glmnet glmnet cv.glmnet
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
#' @param nCores description
#' @param mix description
#' @param filtering_data description
#' @param top_sigs_frac description
#' @param external_essential_genes description
#' @param add_essential_genes description
#' @param return_analysis description
#' @param min_filt_ds description
#' @param predict_res description
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, mix = NULL, ref_type, filtering_data = NULL, lineage_file = NULL, top_genes_frac = 1, medianGEP = TRUE, seed = 123, probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4),
                        sim_fracs = c(seq(0, 0.05, 0.001), seq(0.06, 0.1, 0.005), seq(0.11, 0.25, 0.01)), diff_vals = round(c(log2(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3),
                        min_genes = 3, max_genes = 150, return_sigs = FALSE, return_sigs_filt = FALSE, sigsFile = NULL, minPBcells = 30, minPBsamples = 10, min_filt_ds = 2, predict_res = TRUE,
                        ct_sims = 10,  nCores = 1, top_sigs_frac = 0.05, external_essential_genes = NULL, return_analysis = FALSE, add_essential_genes = TRUE){


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
  out <- normRefMix(ref, mix, ref_type)
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

    if (return_sigs) {
      return(signatures)
    }

    if (!is.null(filtering_data)) {
      message("Filtering signatures by external datasets...")
      filtDS_scored <- scoreFiltDS(ref, labels, filtering_data, ds_cor_cutoff = 0.5, signatures, ncores = nCores)

      if (length(filtDS_scored) > 0) {
        out <- filterSignatures(ref, filtDS_scored, signatures, top_sigs_frac, add_essential_genes, ncores = nCores)
        # out <- filterSignatures2(ref, labels, filtering_data, signatures, top_sigs_frac, add_essential_genes, ncores = nCores)


        if (return_sigs_filt) {
          return(list(all_sigs = signatures,
                      filt_sigs = out$filt_sigs,
                      filt_cts = out$filt_cts,
                      genes_used = rownames(ref)))
        }

        signatures <- out$filt_sigs

      }
    }


    if (!is.null(external_essential_genes)) {
      # TODO: Make a function that add external essential genes
    }


  }else{
    # Load signatures
    message("Loading signatures...")
    signatures <- sigsFile
  }


  # Make simulations
  message("Generating simulations...")
  simulations <- makeSimulations(ref, mix, signatures, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims = ct_sims, ncores = nCores, noise_level = 0.1, seed2use = seed)
  simulations_scored <- scoreSimulations(signatures, simulations, n_sims, sim_fracs, nCores)

  # Learn linear transformation parameters
  # message("Learning linear transformation parameters...")
  # params <- learnParams(simulations_scored, ncores = nCores)
  # data_transfomed <- linearTransform(params, simulations_scored, filtering_data, signatures, min_filt_ds, ref, nCores)

  # Train linear models
  message("Training models...")
  filtDS_scored_filt <- scoreFiltDS(ref, labels, filtering_data, ds_cor_cutoff = 0.5, signatures, ncores = nCores)
  models <- trainModels(simulations_scored, filtering_data, signatures, L1, L2, ncores = nCores)


  # Get spillover matrix
  # message("Generating spillover matrix...")
  # frac2use <- sim_fracs[which.min(abs(sim_fracs - 0.25))]
  # spill_mat <- getSpillOverMat(simulations, signatures_filt, dep_list, models, frac2use)


  # Save results in S4 object
  xCell2Sigs.S4 <- new("xCell2Signatures",
                       signatures = signatures,
                       dependencies = dep_list,
                       params = params,
                       models = models,
                       spill_mat = matrix(),
                       genes_used = rownames(ref))
  message("Custom xCell2.0 reference ready!")


  if (return_analysis) {
    message("Running xCell2Analysis...")
    res <-  xCell2::xCell2Analysis(mix, xcell2sigs = xCell2Sigs.S4, predict = predict_res, spillover = FALSE, ncores = nCores)
    return(res)
  }else{
    return(xCell2Sigs.S4)
  }


}
