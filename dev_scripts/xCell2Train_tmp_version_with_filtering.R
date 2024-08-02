library(tidyverse)

# Load reference
ref.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/bp_ref.rds")
ref = ref.in$ref
labels = ref.in$labels

# Load mixture
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
val_dataset = "GSE20300"
val_tissue = "blood"
mix <- cyto.vals$mixtures[[val_tissue]][[val_dataset]]

# # For mice:
# load("/bigdata/almogangel/xCell2_data/mouse/dataset_petitprez.rda")
# mix <- dataset_petitprez$expr_mat
# human2mouse = TRUE

# Set data type
ref_type = "rnaseq"
val_type = "array"

if (val_type == "array") {
  filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/array_filtering_data.rds")
}else{
  filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/rnaseq_filtering_data.rds")
}

# Remove current validation from filtering data
filt_datasets <- gsub(x = names(filtering_data$mixture), pattern = "#.*", replacement = "")
filtering_data$mixture <- filtering_data$mixture[filt_datasets != val_dataset]
filtering_data$truth <- filtering_data$truth[filt_datasets != val_dataset]


# Load parameters
lineage_file = ref.in$lineage_file
human2mouse = FALSE
seed = 123
num_threads = 40

# For tuning
min_pb_cells = 30
min_pb_samples = 10
min_sc_genes = 1e4
base_count = 3
use_ontology = TRUE
probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.333, 0.4)
diff_vals = round(c(log(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3)
min_genes = 3
max_genes = 200
min_frac_ct_passed = 0.5
return_signatures = FALSE


filter_sigs = FALSE
use_sim2filter = TRUE
sim_fracs = c(0, seq(0.01, 0.25, 0.01))
n_sims = 10
noise_level = 0.2
max_rho_cutoff = 0.5
top_frac_sigs_ds = 0.25
min_ds_frac = 0.5
min_top_genes_frac = 0.5
essen_gene_cutoff = 0.5
min_filt_sigs = 10
top_sigs_frac = 0.05
add_essential_genes = TRUE
return_analysis = FALSE
use_sillover = TRUE



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
sc2pseudoBulk <- function(ref, labels, min_pb_cells, min_pb_samples){


  celltypes <- unique(labels$label)

  groups_list <- lapply(celltypes, function(ctoi){

    ctoi_samples <- labels[labels$label == ctoi,]$sample

    # Calculate maximum possible number of groups given min_pb_cells
    num_groups <- ceiling(length(ctoi_samples) / min_pb_cells)
    if (num_groups < min_pb_samples) {
      num_groups <- min_pb_samples
    }

    # Generate min_pb_samples pseudo samples of CTOI
    if (length(ctoi_samples) > min_pb_samples) {

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
prepRefMix <- function(ref, mix, ref_type, min_sc_genes, base_count, human2mouse){

  if (human2mouse) {
    message("Converting reference genes from human to mouse...")
    data(human_mouse_gene_symbols)
    human_genes <- intersect(rownames(ref), human_mouse_gene_symbols$human)
    ref <- ref[human_genes,]
    rownames(ref) <- human_mouse_gene_symbols[human_genes,]$mouse
  }

  if (ref_type == "sc") {
    message("> Normalizing pseudo bulk reference to CPM.")
    # TODO: For non 10X also do TPM
    lib_sizes <- Matrix::colSums(ref)
    norm_factor <- 1e6 / lib_sizes
    ref_norm <- ref %*% Matrix::Diagonal(x = norm_factor)
    colnames(ref_norm) <- colnames(ref)
    ref <- as.matrix(ref_norm)
    message("> Filtering pseudo bulk genes by variance.")
    genes_var <- sort(apply(ref, 1, var), decreasing = TRUE)
    var_cutoff <- c(1.5, 1, 0.8, 0.5, 0.3, 0.1, 0)
    for (co in var_cutoff) {
      if (co == 0) {
        varGenes2use <- names(genes_var)[1:round(length(genes_var)*0.5)]
      }else{
        varGenes2use <- names(genes_var[genes_var >= co])
      }
      if(length(varGenes2use) > min_sc_genes){
        break
      }else{
        errorCondition("Not enough variable genes in scRNA-Seq reference!")
      }
    }
    ref <- ref[varGenes2use,]

    if (min(ref) < base_count) {
      # Adding base_count to reference restrict inclusion of small changes
      ref <- ref + base_count
    }

  }else{
    if(max(ref) < 50){

      ref <- 2^ref
      if (min(ref) == 1) {
        ref <- ref-1
      }

      if (min(ref) < base_count) {
        # Adding base_count to reference restrict inclusion of small changes
        ref <- ref + base_count
      }

    }else{
      if (min(ref) < base_count) {
        # Adding base_count to reference restrict inclusion of small changes
        ref <- ref + base_count
      }
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
makeQuantiles <- function(ref, labels, probs, num_threads){

  param <- BiocParallel::MulticoreParam(workers = num_threads)
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
createSignatures <- function(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, min_frac_ct_passed, num_threads){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, min_frac_ct_passed){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Set signature cutoffs grid
    param.df <- expand.grid("diff_vals" = diff_vals, "probs" = probs)

    # Find top genes
    type_sigs <- list()

    while (length(type_sigs) < 3 & min_frac_ct_passed >= 0) {

      max_genes_problem <- c()
      for (i in 1:nrow(param.df)){

        # Get a Boolean matrices with genes that pass the quantiles criteria
        diff <- param.df[i, ]$diff_vals # difference threshold
        #lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile cutoff
        lower_prob <- which(as.character(as.numeric(gsub("%", "", rownames(quantiles_matrix[[1]])))/100) == as.character(param.df[i, ]$probs))

        # Sort upper prob gene value for each not_dep_celltypes
        # upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile cutoff
        upper_prob <- which(as.character(as.numeric(gsub("%", "", rownames(quantiles_matrix[[1]])))/100) == as.character(1-param.df[i, ]$probs))
        upper_prob.mat <- sapply(not_dep_celltypes, function(x){
          get(x, quantiles_matrix)[upper_prob,]
        })

        #  Check diff-prob criteria
        diff_genes.mat <- apply(upper_prob.mat, 2, function(x){
          get(type, quantiles_matrix)[lower_prob,] > x + diff
        })

        genes_scores <- apply(diff_genes.mat, 1, function(x){
          names(which(x))
        })
        genes_scores <- genes_scores[lengths(genes_scores) > 0]
        n_ct_passed <- sort(unique(lengths(genes_scores)), decreasing = TRUE)


        # Make signatures
        for (j in n_ct_passed) {

          frac_ct_passed <- round(j/length(not_dep_celltypes), 2)
          if (frac_ct_passed < min_frac_ct_passed) {
            break
          }

          sig_genes <- names(which(lengths(genes_scores) >= j))
          n_genes <- length(sig_genes)

          if (n_genes < min_genes) {
            next
          }

          if (n_genes > max_genes) {
            max_genes_problem <- c(max_genes_problem, TRUE)
            break
          }

          sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, frac_ct_passed, sep = "_")
          type_sigs[[sig_name]] <- sig_genes
        }

      }

      # Fix for cell types that have many DE genes
      if (all(max_genes_problem)) {

        diff_vals_strict <- c(diff_vals, log2(2^max(diff_vals)*2), log2(2^max(diff_vals)*4), log2(2^max(diff_vals)*8), log2(2^max(diff_vals)*16), log2(2^max(diff_vals)*32))
        param.df <- expand.grid("diff_vals" = diff_vals_strict, "probs" = probs)

        for (i in 1:nrow(param.df)){

          # Get a Boolean matrices with genes that pass the quantiles criteria
          diff <- param.df[i, ]$diff_vals # difference threshold
          #lower_prob <- which(probs == param.df[i, ]$probs) # lower quantile cutoff
          lower_prob <- which(as.character(as.numeric(gsub("%", "", rownames(quantiles_matrix[[1]])))/100) == as.character(param.df[i, ]$probs))

          # Sort upper prob gene value for each not_dep_celltypes
          # upper_prob <- nrow(quantiles_matrix[[1]])-lower_prob+1 # upper quantile cutoff
          upper_prob <- which(as.character(as.numeric(gsub("%", "", rownames(quantiles_matrix[[1]])))/100) == as.character(1-param.df[i, ]$probs))
          upper_prob.mat <- sapply(not_dep_celltypes, function(x){
            get(x, quantiles_matrix)[upper_prob,]
          })

          #  Check diff-prob criteria
          diff_genes.mat <- apply(upper_prob.mat, 2, function(x){
            get(type, quantiles_matrix)[lower_prob,] > x + diff
          })

          genes_scores <- apply(diff_genes.mat, 1, function(x){
            names(which(x))
          })
          genes_scores <- genes_scores[lengths(genes_scores) > 0]
          n_ct_passed <- sort(unique(lengths(genes_scores)), decreasing = TRUE)


          # Make signatures
          for (j in n_ct_passed) {

            frac_ct_passed <- round(j/length(not_dep_celltypes), 2)
            if (frac_ct_passed < min_frac_ct_passed) {
              break
            }

            sig_genes <- names(which(lengths(genes_scores) >= j))
            n_genes <- length(sig_genes)

            if (n_genes < min_genes) {
              next
            }

            if (n_genes > max_genes) {
              max_genes_problem <- c(max_genes_problem, TRUE)
              break
            }

            sig_name <-  paste(paste0(type, "#"), param.df[i, ]$probs, diff, n_genes, frac_ct_passed, sep = "_")
            type_sigs[[sig_name]] <- sig_genes
          }

        }

      }



      # Remove duplicate signatures
      type_sigs_sorted <- lapply(type_sigs, function(x) sort(x))
      type_sigs_sorted_collapsed <- sapply(type_sigs_sorted, paste, collapse = ",")
      duplicated_sigs <- duplicated(type_sigs_sorted_collapsed)
      type_sigs <- type_sigs[!duplicated_sigs]

      # Relax parameter until there are at least 3 signatures
      min_frac_ct_passed <- min_frac_ct_passed - 0.05
    }

    return(type_sigs)
  }

  param <- BiocParallel::MulticoreParam(workers = num_threads)
  celltypes <- unique(labels[,2])

  all_sigs <- BiocParallel::bplapply(celltypes, function(type){

    type.sigs <- getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals,
                         min_genes, max_genes, min_frac_ct_passed)


    # Check for minimum 3 signatures per cell type
    if (length(type.sigs) < 3) {
      # Relax prob-diff-min_genes parameters
      probs2use <- c(probs, max(probs)*1.25, max(probs)*1.5, max(probs)*2)
      probs2use <- probs2use[probs2use < 0.5]
      diff_vals2use <- c(min(diff_vals)*0, min(diff_vals)*0.5, min(diff_vals)*0.75, diff_vals)
      min_genes2use <- round(min_genes*0.5)
      min_genes2use <- ifelse(min_genes2use < 5, 5, min_genes2use)

      type.sigs <- getSigs(celltypes, type, dep_list, quantiles_matrix, probs = probs2use, cor_mat, diff_vals = diff_vals2use,
                           min_genes = min_genes2use, max_genes, min_frac_ct_passed = min_frac_ct_passed2use)
    }

    return(type.sigs)
  }, BPPARAM = param)


  all_sigs <- unlist(all_sigs, recursive = FALSE)


  if (length(all_sigs) == 0) {
    stop("No signatures found for reference!")
  }


  return(all_sigs)
}
makeSimulations <- function(ref, mix, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims, noise_level, num_threads){

  getControls <- function(ctoi, controls, gep_mat_linear, sim_fracs, cor_mat, n_sims){

    sampleControls <- function(ctoi, controls, cor_mat){

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
      sampleControls(ctoi, controls, cor_mat)
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


  param <- BiocParallel::MulticoreParam(workers = num_threads)

  celltypes <- unique(labels$label)
  gep_mat_linear <- 2^gep_mat

  sim_list <- lapply(celltypes, function(ctoi){

    # Generate CTOI fractions matrix
    ctoi_mat <- matrix(rep(gep_mat_linear[,ctoi], length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs), dimnames = list(rownames(gep_mat_linear), sim_fracs))
    ctoi_mat_frac <- ctoi_mat %*% diag(sim_fracs)

    # Get control cell types
    if (!is.null(dep_list)) {
      dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
      controls <- celltypes[!celltypes %in% c(ctoi, dep_cts)]
    }else{
      controls <- celltypes[!celltypes == ctoi]
    }

    if (length(controls) == 0) {
      controls <- names(sort(cor_mat[ctoi,])[1])
    }

    # get set of controls that resemble the mixture's shift value
    controls_sets <- getControls(ctoi, controls, gep_mat_linear, sim_fracs, cor_mat, n_sims)

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
          eps_sigma <- mean(apply(mix, 2, sd)) * noise_level
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

  })
  names(sim_list) <- celltypes

  return(sim_list)

}
filterSignatures <- function(shared_cts, shared_genes, labels, filtering_data, simulations, signatures,
                             max_rho_cutoff, top_frac_sigs_ds, min_ds_frac, top_sigs_frac, min_top_genes_frac, essen_gene_cutoff, min_filt_sigs,
                             add_essential_genes, human2mouse, num_threads){


  param <- BiocParallel::MulticoreParam(workers = num_threads)


  celltypes <- unique(labels$label)


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
      if (!is.null(simulations)) {
        filtering_data2use <- simulations[[ctoi]]
      }else{
        return(list(best_sigs = NA,
                    essential_genes = NA))
      }
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
      filter(max_rho >= max_rho_cutoff) %>%
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
      top_frac(top_frac_sigs_ds, wt=rho) %>% # Top 25% correlation per dataset
      #filter(rho >= 0.3) %>%
      group_by(sig) %>%
      summarise(n_sigs = n()) %>%
      mutate(ds_frac = n_sigs/length(ds2use)) %>%
      filter(ds_frac >= min_ds_frac) %>% # Must be in at least 50% of the datasets %>%
      pull(sig) %>%
      unique()

    if (length(top_sigs) == 0) {
      return(list(best_sigs = NA,
                  essential_genes = NA))
    }

    top_genes <- sort(table(unlist(signatures_ctoi[top_sigs])), decreasing = T)/length(top_sigs)
    top_genes <- names(top_genes[top_genes>=min_top_genes_frac]) # Must be in at least 50% of the signatures
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
        top_frac(essen_gene_cutoff, wt=rho_weigted) %>%
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

    top_sigs_frac_adjusted <- ifelse(nrow(rho_weighted_sigs)*top_sigs_frac > min_filt_sigs, top_sigs_frac, min_filt_sigs/nrow(rho_weighted_sigs)) # Minimum 10 signatures

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
learnParams <- function(gep_mat, cor_mat, signatures, dep_list, ref_type, sim_fracs, frac2use, top_spill_value, sc_spill_relaxing_factor, num_threads){

  param <- BiocParallel::MulticoreParam(workers = num_threads)

  celltypes <- colnames(gep_mat)
  gep_mat_linear <- 2^gep_mat


  # Generate mixtures
  mix_list <- BiocParallel::bplapply(celltypes, function(ctoi){

    # Generate CTOI mixture
    ctoi_mat <- matrix(rep(gep_mat_linear[,ctoi], length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs), dimnames = list(rownames(gep_mat_linear), sim_fracs))
    ctoi_mat_frac <- ctoi_mat %*% diag(sim_fracs)

    # Generate control mixture
    if (!is.null(dep_list)) {
      dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
      controls <- celltypes[!celltypes %in% dep_cts]
    }else{
      controls <- celltypes[!celltypes != ctoi]
    }

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
    lm_fit <- lm(sim_fracs ~ scores_transformed)
    m = coef(lm_fit)[[2]]
    n = coef(lm_fit)[[1]]
    # round(scores_transformed*m + n, 2)

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
    if (!is.null(dep_list)) {
      dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
      controls <- celltypes[!celltypes %in% dep_cts]
    }else{
      controls <- celltypes[!celltypes != ctoi]
    }

    controls <- sapply(colnames(cts_mat_frac), function(ct){
      names(sort(cor_mat[controls, ct])[1])
    })
    controls_mat_frac <- gep_mat_linear[,controls] * (1-frac2use)

    # Combine
    mixture <- cts_mat_frac + controls_mat_frac


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
    controls_cts_mat_scores[controls_cts_mat_scores<0] <- 0
    names(controls_cts_mat_scores) <- controls

    final_scores <- round(mix_cts_mat_scores - controls_cts_mat_scores, 2)
    dep_cts <- dep_cts[dep_cts != ctoi]
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
  # TODO: Check why the diagonal is not 0.25
  spill_mat <- spill_mat / diag(spill_mat)


  spill_mat[is.nan(spill_mat)] <- 0
  spill_mat[spill_mat > 1] <- 1

  # TODO: Check this parameter
  if (ref_type == "sc") {
    top_spill_value <- top_spill_value * sc_spill_relaxing_factor
  }

  spill_mat[spill_mat > top_spill_value] <- top_spill_value
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
#' @param min_pb_cells For scRNA-Seq reference only - minimum number of cells in the pseudo-bulk (optional).
#' @param min_pb_samples For scRNA-Seq reference only - minimum number of pseudo-bulk samples (optional).
#' @param min_sc_genes description
#' @param base_count description
#' @param use_ontology A Boolean for using ontological integration (TRUE)
#' @param lineage_file Path to the cell type lineage file generated with `xCell2GetLineage` function (optional).
#' @param probs A numeric vector of probability thresholds to be used for generating signatures (optional).
#' @param diff_vals A numeric vector of delta values to be used for generating signatures (optional).
#' @param min_frac_ct_passed Use for calibration of signatures generation (remove!)
#' @param sim_fracs A vector of mixture fractions to be used in signature filtering (optional).
#' @param min_genes The minimum number of genes to include in the signature (optional).
#' @param max_genes The maximum number of genes to include in the signature (optional).
#' @param sigsFile description
#' @param top_sigs_frac description
#' @param filter_sigs description
#' @param n_sims description
#' @param num_threads description
#' @param human2mouse description
#' @param mix description
#' @param noise_level description
#' @param filtering_data description
#' @param max_rho_cutoff description
#' @param top_frac_sigs_ds description
#' @param min_ds_frac description
#' @param min_top_genes_frac description
#' @param essen_gene_cutoff description
#' @param min_filt_sigs description
#' @param add_essential_genes description
#' @param use_sim2filter description
#' @param top_spill_value description
#' @param sc_spill_relaxing_factor description
#' @param return_analysis description
#' @param return_signatures description
#' @param use_sillover description
#' @param spillover_alpha description
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref,
                        mix = NULL,
                        labels,
                        ref_type,
                        human2mouse = FALSE,
                        lineage_file = NULL,
                        filtering_data = NULL,
                        seed = 123,
                        num_threads = 1,
                        return_signatures = FALSE,
                        return_analysis = FALSE,
                        use_sillover = TRUE,
                        spillover_alpha = 0.2,
                        # For tuning
                        min_pb_cells = 30,
                        min_pb_samples = 10,
                        min_sc_genes = 1e4,
                        base_count = 3,
                        use_ontology = TRUE,
                        probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.333, 0.4),
                        diff_vals = round(c(log(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3),
                        min_genes = 3,
                        max_genes = 200,
                        min_frac_ct_passed = 0.5,
                        top_frac_sigs_ds = 0.25,
                        min_ds_frac = 0.5,
                        min_top_genes_frac = 0.5,
                        essen_gene_cutoff = 0.5,
                        filter_sigs = TRUE,
                        sim_fracs = c(0, seq(0.01, 0.25, 0.01)),
                        n_sims = 10,
                        noise_level = 0.2,
                        max_rho_cutoff = 0.5,
                        min_filt_sigs = 10,
                        top_sigs_frac = 0.05,
                        add_essential_genes = TRUE,
                        use_sim2filter = TRUE,
                        top_spill_value = 0.5,
                        sc_spill_relaxing_factor = 0.5
){

  set.seed(seed)

  # Validate inputs
  inputs_validated <- validateInputs(ref, labels, ref_type)
  ref <- inputs_validated$ref
  labels <- inputs_validated$labels

  # Generate pseudo bulk from scRNA-Seq reference
  if (ref_type == "sc") {
    message("Converting scRNA-seq reference to pseudo bulk...")
    pb_data <- sc2pseudoBulk(ref, labels, min_pb_cells, min_pb_samples)
    ref <- pb_data$ref
    labels <- pb_data$labels
  }

  # Prepare reference and mixture data: human to mouse genes transformation, normalization, base_counts, log2 transformation, shared genes
  out <- prepRefMix(ref, mix, ref_type, min_sc_genes, base_count, human2mouse)
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
  quantiles_matrix <- makeQuantiles(ref, labels, probs, num_threads)
  message("Generating signatures...")
  signatures <- createSignatures(labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, min_frac_ct_passed, num_threads)

  if (return_signatures) {
    xCell2.S4 <- new("xCell2Object",
                     signatures = signatures,
                     dependencies = list(),
                     params = data.frame(),
                     spill_mat = matrix(),
                     genes_used = shared_genes)
    return(xCell2.S4)
  }

  if (filter_sigs) {
    # Generate simulations
    message("Generating simulations...")
    if (use_sim2filter) {
      simulations <- makeSimulations(ref, mix, labels, gep_mat, ref_type, dep_list, cor_mat, sim_fracs, n_sims, noise_level, num_threads)
    }else{
      simulations <- NULL
    }

    # Filter signatures
    message("Filtering signatures...")
    # Identify cell types that exist in the filtering datasets
    shared_cts <- intersect(unique(labels$label), unlist(sapply(filtering_data$truth, rownames)))
    signatures <- filterSignatures(shared_cts, shared_genes, labels, filtering_data, simulations, signatures,
                                   max_rho_cutoff, top_frac_sigs_ds, min_ds_frac, top_sigs_frac, min_top_genes_frac, essen_gene_cutoff, min_filt_sigs,
                                   add_essential_genes, human2mouse, num_threads)
  }

  # Learn linear transformation parameters
  message("Learning linear transformation and spillover parameters...")
  params <- learnParams(gep_mat, cor_mat, signatures, dep_list, ref_type, sim_fracs, frac2use = 0.25, top_spill_value, sc_spill_relaxing_factor, num_threads)
  spill_mat <- params$spillmat
  params <- params$params


  # Save results in S4 object
  xCell2.S4 <- new("xCell2Object",
                   signatures = signatures,
                   dependencies = dep_list,
                   params = params,
                   spill_mat = spill_mat,
                   genes_used = shared_genes)

  message("Your custom xCell2 reference object is ready!")
  message("> Please consider sharing your xCell2 reference with others here: https://dviraran.github.io/xCell2ref")


  if (return_analysis) {
    message("Running xCell2Analysis...")
    res <- xCell2::xCell2Analysis(mix, xcell2object = xCell2.S4, spillover = use_sillover, spillover_alpha = spillover_alpha, num_threads = num_threads)
    return(res)
  }else{
    return(xCell2.S4)
  }

}
