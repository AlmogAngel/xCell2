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
prepRefMix <- function(ref, mix, ref_type, human2mouse){

  if (human2mouse) {
    message("Converting reference genes from human to mouse...")

    data(human_mouse_gene_symbols)

    human_genes <- intersect(rownames(ref), human_mouse_gene_symbols$human)
    ref <- ref[human_genes,]
    rownames(ref) <- human_mouse_gene_symbols[human_genes,]$mouse

  }


  if (ref_type == "sc") {
    message("Normalizing pseudo bulk reference to CPM.")
    # TODO: For non 10X also do TPM
    lib_sizes <- Matrix::colSums(ref)
    norm_factor <- 1e6 / lib_sizes
    ref_norm <- ref %*% Matrix::Diagonal(x = norm_factor)
    colnames(ref_norm) <- colnames(ref)
    ref <- as.matrix(ref_norm) # TODO: Find a way to reduce memory usage by keeping matrix sparse
    ref <- ref + 1
  }else{
    # Adding 3 to reference restrict inclusion of small changes
    if(max(ref) < 50 & !is.null(mix)){
      # Unlog2
      ref <- 2^ref
      if (min(ref) == 1) {
        ref <- ref-1
      }
      ref <- ref + 3
    }else{
      ref <- ref + 3
    }
  }

  # log2-transformation
  ref <- log2(ref)


  if(max(mix) >= 50 & !is.null(mix)){
    mix <- log2(mix+1)
  }

  if (!is.null(mix)) {
    # Select shared genes
    shared_genes <- intersect(rownames(ref), rownames(mix))
    message(paste0(length(shared_genes), " genes are shared between reference and mixture."))
    ref <- ref[shared_genes,]
    mix <- mix[shared_genes,]
  }


  return(list(ref.out = ref, mix.out = mix))
}
makeGEPMat <- function(ref, labels){

  celltypes <- unique(labels$label)

  gep_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      type_vec <- Rfast::rowMedians(as.matrix(ref[,type_samples]))
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
    param.df <- param.df[param.df$diff_vals != 0 | param.df$probs <= 0.05,] # Remove diff-prob combination that is too permissive...
    ntop <- round(length(not_dep_celltypes)*top_genes_frac) # Number of top values given top_genes_frac of the not dependent cell types
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
makeSimulations <- function(ref, mix, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims, ncores, noise_level, seed2use){

  getControls <- function(ctoi, controls, gep_mat_linear, dep_cts, sim_fracs, cor_mat, n_sims){

    sampleControls <- function(ctoi, controls, dep_cts, cor_mat){

      cts_cors <- sort(cor_mat[ctoi, controls], decreasing = TRUE)
      controls <- controls[!controls %in% names(which(cts_cors > 0.8))]

      # Sample 3-9 non dependent control cell types (if possible)
      if (length(controls) > 3) {
        controls2use <- sample(controls, sample(3:min(8, length(controls)), 1), replace = FALSE)
      }else{
        controls2use <- controls
      }

      # # For half of the simulations add a related cell type(s) to mimic some spillover effect
      # related_cts <- names(cts_cors[1:round(length(cts_cors)*0.2)]) # Top 20% related cell type
      #
      # if (length(related_cts) > 0) {
      #   if (runif(1) <= 0.5) {
      #     if (length(related_cts) < 4) {
      #       controls2use <- unique(c(controls2use, sample(related_cts, 1)))
      #
      #     }else{
      #       controls2use <- unique(c(controls2use, sample(related_cts, 2)))
      #     }
      #   }
      # }
      return(controls2use)

    }

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
    controls_sets <- getControls(ctoi, controls, gep_mat_linear, dep_cts, sim_fracs, cor_mat, n_sims)

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

        if (is.null(mix)) {
          eps_sigma <- 2 * noise_level
        }else{
          eps_sigma <- sd(as.vector(mix)) * noise_level
        }
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
filterSignatures <- function(shared_genes, labels, filtering_data, simulations, signatures, top_sigs_frac, add_essential_genes, human2mouse, ncores){


  param <- BiocParallel::MulticoreParam(workers = ncores)


  celltypes <- unique(labels$label)
  shared_cts <- intersect(celltypes, unlist(sapply(filtering_data$truth, rownames)))


  if (human2mouse) {
    data(human_mouse_gene_symbols)

    mouse_genes_mixtures <- lapply(filtering_data$mixture, function(m){
      human_genes <- intersect(rownames(m), human_mouse_gene_symbols$human)
      m <- m[human_genes,]
      rownames(m) <- human_mouse_gene_symbols[human_genes,]$mouse
      return(m)
    })
    filtering_data$mixture <- mouse_genes_mixtures
  }


  filt_sigs <- BiocParallel::bplapply(celltypes, function(ctoi){

    # Get CTOI signatures
    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

    # Use external filtering data/simulations
    if (ctoi %in% shared_cts) {
      ds2use <- names(which(sapply(filtering_data$truth, function(x){ctoi %in% rownames(x)})))
      filtering_data2use <- filtering_data$mixture[ds2use]
    }else{
      filtering_data2use <- simulations[[ctoi]]
    }

    # Calculate CTOI correlations for each dataset in the filtering data
    ds_cors_list <- lapply(names(filtering_data2use), function(ds){

      ctoi_mix <- filtering_data2use[[ds]]

      if (ctoi %in% shared_cts) {
        ctoi_samples <- filtering_data$truth[[ds]][ctoi,]
        ctoi_samples <- names(ctoi_samples[!is.na(ctoi_samples)])
        ctoi_samples <- ctoi_samples[ctoi_samples %in% colnames(ctoi_mix)]
        genes2use <- intersect(rownames(ctoi_mix), shared_genes)
        ctoi_mix <- ctoi_mix[genes2use, ctoi_samples]
      }else{
        genes2use <- shared_genes
      }


      #  Rank filtering dataset
      ctoi_filt_ds_ranked <- singscore::rankGenes(ctoi_mix)

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
      if (ctoi %in% shared_cts) {
        fracs <- filtering_data$truth[[ds]][ctoi, ctoi_samples]
      }else{
        fracs <- as.numeric(gsub(".*%%", "", colnames(ctoi_mix)))
      }

      cors <- apply(scores, 2, function(x){
        cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
      })

      return(cors)

    })
    names(ds_cors_list) <- names(filtering_data2use)

    # Get number of samples per dataset
    if (ctoi %in% shared_cts) {
      ds2n_samples <- enframe(sapply(names(ds_cors_list), function(ds){
        mix_samples <- colnames(filtering_data$mixture[[ds]])
        truth_samples <- filtering_data$truth[[ds]][ctoi,]
        truth_samples <- names(truth_samples[!is.na(truth_samples)])
        n_samples <- length(intersect(mix_samples, truth_samples))
        n_samples
      }), name = "ds", value = "n_samples")
    }else{
      ds2n_samples <- enframe(sapply(names(ds_cors_list), function(ds){
        ncol(filtering_data2use[[ds]])
      }), name = "ds", value = "n_samples")
    }


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
    top_genes <- intersect(top_genes, shared_genes)

    if (length(top_genes) > 0) {
      ds_top_genes_cors <- lapply(ds2use, function(ds){

        if (ctoi %in% shared_cts) {
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

        }else{
          top_genes_ctoi_mat <- simulations[[ctoi]][[ds]][top_genes,]
          true_fracs <- as.numeric(gsub(".*%%", "", colnames(top_genes_ctoi_mat)))
        }



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
  names(filt_sigs) <- celltypes


  # Filter genes
  for(ctoi in celltypes){
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

  # out <- list(filt_sigs = signatures,
  #             filt_cts = filts_ct)

  return(signatures)
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
    # ctoi_sims <- Reduce(rbind, ctoi_sims)
    names(ctoi_sims) <- paste0("sim-", 1:length(ctoi_sim_list))

    ctoi_sims
  }, BPPARAM = param)
  names(sims_scored) <- celltypes


  return(sims_scored)

}
learnParams <- function(simulations_scored, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)

  celltypes <- names(simulations_scored)

  linearParams <- BiocParallel::bplapply(celltypes, function(ctoi){


    tps <- Rfast::rowMedians(sapply(simulations_scored[[ctoi]], function(sim){

      frac <- sim[,ncol(sim)]
      mean_scores <- rowMeans(sim[,-ncol(sim)])

      tp <- try(minpack.lm::nlsLM(mean_scores ~ a * frac^b, start = list(a=1, b=1), control = list(maxiter = 500)), silent = TRUE)
      a = coef(tp)[[1]]
      b = coef(tp)[[2]]

      mean_scores_transformed <- (mean_scores^(1/b)) / a
      lm_fit <- lm(frac ~ mean_scores_transformed)
      m = coef(lm_fit)[[2]]
      n = coef(lm_fit)[[1]]


      return(c(a = a, b = b, m = m, n = n))

    }))


    return(tibble(celltype = ctoi, a = tps[[1]], b = tps[[2]], m = tps[[3]], n = tps[[4]]))

  }, BPPARAM = param)


  bind_rows(linearParams) %>%
    return(.)

}
getSpillOverMat <- function(gep_mat, cor_mat, signatures, dep_list, params, frac2use, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)


  celltypes <- colnames(gep_mat)
  gep_mat_linear <- 2^gep_mat
  if (round(min(gep_mat_linear)) == 3) {
    gep_mat_linear <- gep_mat_linear-3
  }


  # Make ctoi and control mixtures
  ctoi_mat_frac <- gep_mat_linear * frac2use
  control_mat <- sapply(celltypes, function(ctoi){
    control <- names(sort(cor_mat[,ctoi]))[1]
    gep_mat_linear[,control]
  })
  control_mat_frac <- control_mat*(1-frac2use)
  mix_mat <- ctoi_mat_frac + control_mat_frac

  mix_mat_ranked <- singscore::rankGenes(mix_mat)
  control_mat_ranked <- singscore::rankGenes(control_mat)


  ctoi_spill_scores <- BiocParallel::bplapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]
    dep_cts <- unlist(dep_list[[ctoi]])
    a <- pull(filter(params, celltype == ctoi), a)
    b <- pull(filter(params, celltype == ctoi), b)

    mix_scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    mix_scores <- Rfast::rowmeans(mix_scores)
    mix_scores <- (mix_scores^(1/b)) / a
    names(mix_scores) <- colnames(mix_mat_ranked)


    control_scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(control_mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    control_scores <- Rfast::rowmeans(control_scores)
    control_scores <- (control_scores^(1/b)) / a
    names(control_scores) <- colnames(control_mat_ranked)

    final_scores <- mix_scores - control_scores
    final_scores[dep_cts] <- 0

    return(final_scores)

  }, BPPARAM = param)
  names(ctoi_spill_scores) <- celltypes

  # Rows are ctoi signatures scores
  # Columns are cell type in mixtures
  spill_mat <- Reduce(rbind, ctoi_spill_scores)
  rownames(spill_mat) <- celltypes


  # Clean and normalize spill matrix
  spill_mat[spill_mat < 0] <- 0
  spill_mat <- spill_mat / diag(spill_mat)
  spill_mat[is.nan(spill_mat)] <- 0
  spill_mat[spill_mat > 1] <- 1


  # TODO: Check this parameter
  spill_mat[spill_mat > 0.5] <- 0.5
  diag(spill_mat) <- 1

  # pheatmap::pheatmap(spill_mat, cluster_rows = F, cluster_cols = F)

  return(spill_mat)

}
scoreFiltDS <- function(shared_genes, labels, filtering_data, ds_cor_cutoff, min_ds2use, signatures, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)

  filt_cts <- intersect(unique(labels$label), unlist(sapply(filtering_data$truth, rownames)))
  stableG <- singscore::getStableGenes(5, type = 'blood')


  if (length(filt_cts) == 0) {
    warningCondition("Not cell types found for filtering signatures.")
    return(NA)
  }

  ct_ds_scores_list <- BiocParallel::bplapply(filt_cts, function(ctoi){

    # Get CTOI signatures
    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

    # Score
    ds2use <- names(which(unlist(lapply(filtering_data$truth, function(ds){ctoi %in% rownames(ds)}))))
    filtering_data_ctoi <- list("mixture" = filtering_data$mixture[ds2use],
                                "truth" = filtering_data$truth[ds2use])

    genes2use <- intersect(Reduce(intersect, lapply(filtering_data_ctoi$mixture, rownames)), shared_genes)

    ds_scores_list <- lapply(ds2use, function(ds){

      ctoi_mix <- filtering_data$mixture[[ds]]
      ctoi_samples <- filtering_data$truth[[ds]][ctoi,]
      ctoi_samples <- names(ctoi_samples[!is.na(ctoi_samples)])
      ctoi_samples <- ctoi_samples[ctoi_samples %in% colnames(ctoi_mix)]

      #  Rank filtering dataset
      ctoi_filt_ds_ranked <- singscore::rankGenes(ctoi_mix[genes2use, ctoi_samples], stableGenes = stableG)

      # Get stable genes score
      median(singscore::simpleScore(ctoi_filt_ds_ranked, upSet = stableG, centerScore = FALSE)$TotalScore)

      # Score
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        suppressWarnings(singscore::simpleScore(ctoi_filt_ds_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
      })


      # Remove signatures that failed to score
      if (class(scores)[1] == "list") {
        scores <- scores[lengths(scores) != 0]
        scores <- sapply(scores, cbind)
      }

      scores[scores < 0] <- 0 # Scores might be negative if all signature's genes are zero


      fracs <- filtering_data$truth[[ds]][ctoi, ctoi_samples]

      return(cbind(scores, "frac" = fracs))

    })
    names(ds_scores_list) <- ds2use

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

    if (length(ds2use) < min_ds2use) {
      return(NA)
    }

    ds_scores_list <- ds_scores_list[ds2use]

    return(ds_scores_list)
    #final_training_ds <- Reduce(rbind, ds_scores_list)

  }, BPPARAM = param)
  names(ct_ds_scores_list) <- filt_cts

  ct_ds_scores_list <- ct_ds_scores_list[sapply(ct_ds_scores_list, function(x)!is.na(x)[1])]

  return(ct_ds_scores_list)

}
trainModels <- function(simulations_scored, params, alpha, ncores, seed2use){

  fitModel <- function(data, a, b, alpha, nCores){

    # Each filtering data / simulation get a different model
    models_coef <- lapply(data, function(ds){

      Y <- ds[,ncol(ds)]
      X <- ds[,-ncol(ds)]

      # Linear transformation
      X <- (X^(1/b)) / a

      # # Scale
      # X_means <- apply(X, 2, mean)
      # X_sds <- apply(X, 2, sd)
      # X <- scale(X, center = X_means, scale = X_sds)


      num_samples <- nrow(X)

      # Strategy for determining 'nfold'
      if (num_samples < 30) {
        nfold <- num_samples
        grouped <- FALSE
      } else if (num_samples >= 30 && num_samples <= 100) {
        nfold <- min(20, num_samples)
        grouped <- TRUE
      } else {
        nfold <- 10
        grouped <- TRUE
      }

      cv_fit <- glmnet::cv.glmnet(X, Y, nfolds = nfold, grouped = grouped, alpha = alpha, family = "gaussian")
      models_coef <- as.matrix(coef(cv_fit, s = "lambda.min"))

      intercept <- models_coef[1,]
      betas <- models_coef[-1,]

      # Scale and predict
      new_X <- scale((ds[,-ncol(ds)]^(1/b)) / a , center = X_means, scale = X_sds)
      (new_X %*% betas) + intercept
      Y

      return(tibble(X_means = list(X_means), X_sds = list(X_sds), betas = list(betas), intercept = list(intercept)))
    })
    models_coef <- bind_rows(models_coef)




    return(models_coef)

  }

  set.seed(seed2use)
  param <- BiocParallel::MulticoreParam(workers = ncores)


  celltypes <- names(simulations_scored)

  models_params_list <- BiocParallel::bplapply(celltypes, function(ctoi){

    a <- pull(filter(params, celltype == ctoi), a)
    b <- pull(filter(params, celltype == ctoi), b)
    data <- simulations_scored[[ctoi]]

    model_params <- fitModel(data, a = a, b = b, alpha = alpha, nCores = ncores)

    return(tibble(celltype = ctoi, models_params = list(model_params)))

  }, BPPARAM = param)


  bind_rows(models_params_list) %>%
    left_join(params, by = "celltype") %>%
    return(.)

}


learnParams <- function(gep_mat, cor_mat, signatures, dep_list, params, sim_fracs, frac2use, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)

  celltypes <- colnames(gep_mat)
  gep_mat_linear <- 2^gep_mat
  if (round(min(gep_mat_linear)) == 1) {
    gep_mat_linear <- gep_mat_linear-1
  }

  # Generate mixtures
  mix_list <- BiocParallel::bplapply(celltypes, function(ctoi){

    # Generate CTOI mixture
    ctoi_mat <- matrix(rep(gep_mat_linear[,ctoi], length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs), dimnames = list(rownames(gep_mat_linear), sim_fracs))
    ctoi_mat_frac <- ctoi_mat %*% diag(sim_fracs)

    # Generate control mixture
    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
    controls <- celltypes[!celltypes %in% dep_cts]
    control <- names(sort(cor_mat[controls,ctoi])[1])
    controls_mat <- matrix(rep(gep_mat_linear[,control], length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs), dimnames = list(rownames(gep_mat_linear), sim_fracs))
    controls_mat_frac <- controls_mat %*% diag(1-sim_fracs)

    # Combine
    mixture <- ctoi_mat_frac + controls_mat_frac
    colnames(mixture) <- paste0(ctoi, "^^", control, "%%", sim_fracs)

    return(mixture)

  }, BPPARAM = param)
  names(mix_list) <- celltypes

  # Learn linear parameters
  linear_params <- BiocParallel::bplapply(celltypes, function(ctoi){

    # Get scores
    mix_mat_ranked <- singscore::rankGenes(mix_list[[ctoi]])
    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]
    scores <- rowMeans(sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    }))
    # plot(scores, sim_fracs)

    # Get transformation parameters
    tp <- try(minpack.lm::nlsLM(scores ~ a * sim_fracs^b, start = list(a=1, b=1), control = list(maxiter = 500)), silent = TRUE)
    a = coef(tp)[[1]]
    b = coef(tp)[[2]]
    # plot((scores^(1/b)) / a, sim_fracs)

    # Get linear model parameters
    scores_transformed <- (scores^(1/b)) / a
    lm_fit <- lm(frac ~ scores_transformed)
    m = coef(lm_fit)[[2]]
    n = coef(lm_fit)[[1]]


    return(tibble(celltype = ctoi, a = a, b = b, m = m, n = n))
  }, BPPARAM = param)
  linear_params <- bind_rows(linear_params)

  # Learn spillover parameters
  ctoi_spill_scores <- BiocParallel::bplapply(celltypes, function(ctoi){

    signatures_ctoi <- signatures[gsub("#.*", "", names(signatures)) %in% ctoi]
    a <- pull(filter(linear_params, celltype == ctoi), a)
    b <- pull(filter(linear_params, celltype == ctoi), b)
    m <- pull(filter(linear_params, celltype == ctoi), m)
    n <- pull(filter(linear_params, celltype == ctoi), n)

    # Generate frac2use mixtures
    cts_mat_frac <- gep_mat_linear * frac2use

    # Generate control mixture
    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
    dep_cts <- dep_cts[dep_cts != ctoi]
    controls <- celltypes[!celltypes %in% dep_cts]
    controls <- sapply(colnames(cts_mat_frac), function(ct){
      names(sort(cor_mat[controls, ct])[1])
    })
    controls_mat_frac <- gep_mat_linear[,controls] * (1-frac2use)

    # Combine
    mixture <- ctoi_mat_frac + controls_mat_frac


    # Get results for all cell type mixtures
    mix_cts_mat_ranked <- singscore::rankGenes(mixture)
    mix_cts_mat_scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_cts_mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    mix_cts_mat_scores <- Rfast::rowmeans(mix_cts_mat_scores)
    mix_cts_mat_scores <- (mix_cts_mat_scores^(1/b)) / a
    mix_cts_mat_scores <- mix_cts_mat_scores*m + n
    names(mix_cts_mat_scores) <- colnames(mix_cts_mat_ranked)


    # Get results for all cell type controls
    controls_cts_mat_ranked <- singscore::rankGenes(controls_mat_frac)
    colnames(controls_cts_mat_ranked) <- make.unique(colnames(controls_cts_mat_ranked))
    controls_cts_mat_scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(controls_cts_mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    controls_cts_mat_scores[controls_cts_mat_scores<0] <- 0
    controls_cts_mat_scores <- Rfast::rowmeans(controls_cts_mat_scores)
    controls_cts_mat_scores <- (controls_cts_mat_scores^(1/b)) / a
    controls_cts_mat_scores <- controls_cts_mat_scores*m + n


    final_scores <- round(mix_cts_mat_scores - controls_cts_mat_scores, 4)
    final_scores[dep_cts] <- 0
    final_scores[final_scores < 0] <- 0



    return(final_scores)

  }, BPPARAM = param)
  names(ctoi_spill_scores) <- celltypes

  # Rows are ctoi signatures scores
  # Columns are cell type in mixtures
  spill_mat <- Reduce(rbind, ctoi_spill_scores)
  rownames(spill_mat) <- celltypes

  # Clean and normalize spill matrix
  spill_mat[spill_mat < 0] <- 0
  spill_mat <- spill_mat / diag(spill_mat)
  spill_mat[is.nan(spill_mat)] <- 0
  spill_mat[spill_mat > 1] <- 1

  # TODO: Check this parameter
  spill_mat[spill_mat > 0.5] <- 0.5
  diag(spill_mat) <- 1

  # pheatmap::pheatmap(spill_mat, cluster_rows = F, cluster_cols = F)


  return(list(params = linear_params, spillmat = spill_mat))
}

#' @slot signatures list of xCell2 signatures
#' @slot dependencies list of cell type dependencies
#' @slot params data frame of cell type linear transformation parameters
#' @slot spill_mat matrix of cell types spillover
#' @slot genes_used character vector of genes names used to train the signatures
#' @importFrom methods new
# Create S4 object for the new reference
setClass("xCell2Object", slots = list(
  signatures = "list",
  dependencies = "list",
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
#' @importFrom minpack.lm nlsLM
#' @importFrom glmnet cv.glmnet
#' @importFrom Rfast rowMedians rowmeans rowsums Sort
#' @importFrom Matrix rowMeans rowSums colSums
#' @importFrom singscore rankGenes simpleScore
#' @param ref A reference gene expression matrix (genes in rows samples/cells in columns).
#' @param labels A data frame in which the rows correspond to samples in the ref. The data frame must have four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545" or NA if there is no ontology).
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample that match the column name in ref.
#'   "dataset": sample's source dataset or subject (for single-cell).
#' @param ref_type Gene expression data type: "rnaseq" for bulk RNA-Seq, "array" for micro-array, or "sc" for scRNA-Seq.
#' @param seed Set seed for reproducible results (optional).
#' @param minPBcells For scRNA-Seq reference only - minimum number of cells in the pseudo-bulk (optional).
#' @param minPBgroups For scRNA-Seq reference only - minimum number of pseudo-bulk samples (optional).
#' @param use_ontology A Boolean for using ontological integration (TRUE)
#' @param lineage_file Path to the cell type lineage file generated with `xCell2GetLineage` function (optional).
#' @param probs A numeric vector of probability thresholds to be used for generating signatures (optional).
#' @param diff_vals A numeric vector of delta values to be used for generating signatures (optional).
#' @param top_genes_frac Use for calibration of signatures generation (remove!)
#' @param sim_fracs A vector of mixture fractions to be used in signature filtering (optional).
#' @param min_genes The minimum number of genes to include in the signature (optional).
#' @param max_genes The maximum number of genes to include in the signature (optional).
#' @param sigsFile description
#' @param top_sigs_frac description
#' @param return_sigs description
#' @param return_sigs_filt description
#' @param ct_sims description
#' @param nCores description
#' @param human2mouse description
#' @param mix description
#' @param filtering_data description
#' @param external_essential_genes description
#' @param add_essential_genes description
#' @param return_analysis description
#' @param regAlpha description
#' @param use_sillover description
#' @param predict_res description
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, mix = NULL, ref_type, seed = 123, nCores = 1, human2mouse = FALSE,
                        use_ontology = TRUE, lineage_file = NULL, top_genes_frac = 1,
                        probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4),
                        diff_vals = round(c(log2(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3),
                        min_genes = 3,
                        max_genes = 150,
                        add_essential_genes = TRUE,
                        filtering_data = NULL,
                        sim_fracs = c(seq(0, 0.05, 0.001), seq(0.06, 0.1, 0.005), seq(0.11, 0.25, 0.01)),
                        return_sigs = FALSE, return_sigs_essen = FALSE, return_sigs_filt = FALSE, sigsFile = NULL, minPBcells = 30, minPBsamples = 10, regAlpha = 0.5, predict_res = TRUE,
                        ct_sims = 10, external_essential_genes = NULL, return_analysis = FALSE, use_sillover = TRUE, top_sigs_frac = 0.05){


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
  out <- prepRefMix(ref, mix, ref_type, human2mouse)
  ref <- out$ref.out
  mix <- out$mix.out
  shared_genes <- rownames(ref)

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  gep_mat <- makeGEPMat(ref, labels)
  cor_mat <- getCellTypeCorrelation(gep_mat, ref_type)

  # Get cell type dependencies list
  if (use_ontology) {
    message("Loading dependencies...")
    if (is.null(lineage_file)) {
      dep_list <- xCell2::xCell2GetLineage(labels, out_file = NULL)
    }else{
      dep_list <- getDependencies(lineage_file)
    }
  }else{
    message("Skipping ontological integration")
    dep_list <- NULL
  }

  # Generate signatures
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs, ncores = nCores)
  message("Generating signatures...")
  signatures <- createSignatures(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, top_genes_frac, ncores = nCores)

  # Generate simulations
  message("Generating simulations...")
  simulations <- makeSimulations(ref, mix, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims = ct_sims, ncores = nCores, noise_level = 0.2, seed2use = seed)

  # Filter signatures
  message("Filtering signatures...")
  signatures <- filterSignatures(shared_genes, labels, filtering_data, simulations, signatures, top_sigs_frac, add_essential_genes, human2mouse, ncores = nCores)

  # Score simulations
  #message("Scoring simulations for parameter tuning...")
  #simulations_scored <- scoreSimulations(signatures, simulations, n_sims, sim_fracs, nCores)

  # Learn linear transformation parameters
  message("Learning linear transformation and spillover parameters...")
  params <- learnParams(gep_mat, cor_mat, signatures, dep_list, params, sim_fracs, frac2use = 0.25, nCores)
  spill_mat <- params$spillmat
  params <- params$params

  # # Get spillover parameters
  # if (use_sillover) {
  #   message("Generating spillover matrix...")
  #   spill_mat <- getSpillOverMat(gep_mat, cor_mat, signatures, dep_list, params, frac2use = 0.25, nCores)
  # }else{
  #   spill_mat <- matrix
  # }


  # Learn scores transformation parameters
  # TODO: Parameters that transform signatures scores to estimated cell fractions

  # # Train linear models
  # message("Learning signatures regularization parameters...")
  # filtDS_scored <- scoreFiltDS(shared_genes, labels, filtering_data, ds_cor_cutoff = 0.5, min_ds2use = 1, signatures, ncores = nCores)
  # params <- trainModels(simulations_scored, params, alpha = 0, ncores = nCores, seed2use = seed)



  # Save results in S4 object
  xCell2.S4 <- new("xCell2Object",
                   signatures = signatures,
                   dependencies = dep_list,
                   params = params,
                   spill_mat = spill_mat,
                   genes_used = rownames(ref))

  message("Your custom xCell2 reference object is ready!")
  message("> You can help others by sharing your reference here: https://dviraran.github.io/xCell2ref")


  if (return_analysis) {
    message("Running xCell2Analysis...")
    res <- xCell2::xCell2Analysis(mix, xcell2object = xCell2.S4, spillover = use_sillover, ncores = nCores)
    return(res)
  }else{
    return(xCell2.S4)
  }

}
