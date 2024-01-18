library(tidyverse)


# Load reference
ref.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/kass_blood_ref.rds")
ref = ref.in$ref
labels = ref.in$labels

# Load mixture
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
val_dataset = "GSE65133"
val_tissue = "blood"
mix <- cyto.vals$mixtures[[val_tissue]][[val_dataset]]

# Set data type
ref_type = "rnaseq"

# Load filtering data
filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/blood_array_filtering_data.rds")

# Remove current validation from filtering data
filt_datasets <- gsub(x = colnames(filtering_data$mixture), pattern = "#.*", replacement = "")
filtering_data$mixture <- filtering_data$mixture[,filt_datasets != val_dataset]
filtering_data$truth <- filtering_data$truth[,filt_datasets != val_dataset]
filtering_data$truth <- filtering_data$truth[!apply(filtering_data$truth, 1, function(x){all(is.na(x))}),]

# Load parameters
lineage_file = ref.in$lineage_file
sim_fracs = c(0, seq(0.01, 0.25, 0.01), seq(0.3, 1, 0.05))
diff_vals = round(c(log2(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3)
probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4)
min_genes = 3
max_genes = 150
return_sigs = FALSE
sigsFile = NULL
minPBcells = 30
minPBsamples = 10
top_genes_frac = 1
medianGEP = TRUE
sim_noise = NULL
ct_sims = 10
samples_frac = 0.1
seed = 123
nCores = 20
simMethod = "mix_thin"
# c("ref_multi", "mix_thin")


# Load signature
# sigsFile = "/bigdata/almogangel/xCell2_data/dev_data/sigs/BG_blood_ts_blood_sigs.rds"



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
filterSignatures <- function(ref, labels, filtering_data, signatures, n_sigs, ncores){

  param <- BiocParallel::MulticoreParam(workers = ncores)

  normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }

  celltypes <- unique(labels$label)
  shared_cts <- intersect(celltypes, rownames(filtering_data$truth))

  filt_sigs <- BiocParallel::bplapply(shared_cts, function(ctoi){

    # Get CTOi signatures
    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

    # Use only samples with truth values
    ctoi_samples <- colnames(filtering_data$truth)[!is.na(filtering_data$truth[ctoi,])]

    # Get filtering data
    ctoi_filt_data <- filtering_data$mixture[,ctoi_samples]

    # Calculate CTOI for each dataset in the filtering data
    ctoi_datasets <- gsub("#.*", "", ctoi_samples)
    ds_cors_list <- lapply(unique(ctoi_datasets), function(ds){

      # Get dataset's samples
      ctoi_filt_data_ds <- ctoi_filt_data[,ctoi_datasets == ds]

      #  Rank dataset
      mix_ranked <- singscore::rankGenes(ctoi_filt_data_ds)

      # Score and normalize
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){

        suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
      })
      scores <- apply(scores, 2, function(sample){
        normalize(sample)
      })

      # Calculate correlations
      fracs <- filtering_data$truth[ctoi, colnames(mix_ranked)]
      cors <- apply(scores, 2, function(x){
        cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
      })

      return(sort(cors, decreasing = T))

    })
    names(ds_cors_list) <- unique(ctoi_datasets)

    # Rank signatures by dataset
    ranked_sigs <- lapply(ds_cors_list, function(ds) {
      rank(ds)
    })

    best_sigs <- enframe(ranked_sigs, name = "dataset") %>%
      unnest_longer(value, values_to = "rank", indices_to = "sig") %>%
      group_by(sig) %>%
      summarise(rank = mean(rank)) %>%
      top_n(n_sigs, wt = rank) %>%
      pull(sig)


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
makeSimulations <- function(ref, labels, mix, gep_mat, sim_fracs, ncores, seed2use){

  set.seed(seed2use)
  param <- BiocParallel::MulticoreParam(workers = ncores)

  # Use QN
  y <- cbind(ref, mix)
  y <- limma::normalizeBetweenArrays(y)
  ref.norm <- y[,1:ncol(ref)]

  celltypes <- unique(labels$label)

  ref.norm.linear <- 2^ref.norm
  if (round(min(ref.norm.linear)) == 1) {
    ref.norm.linear <- ref.norm.linear-1
  }

  sim_list <- BiocParallel::bplapply(celltypes, function(ctoi){

    ctoi_vec <- gep_mat[,ctoi]
    ctoi_vec_linear <- 2^ctoi_vec
    if (min(ctoi_vec_linear) == 1) {
      ctoi_vec_linear <- ctoi_vec_linear-1
    }

    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
    controls <- celltypes[!celltypes %in% dep_cts]
    controls_vec <- Rfast::rowMedians(gep_mat[,controls])
    controls_vec_linear <- 2^controls_vec
    if (min(controls_vec_linear) == 1) {
      controls_vec_linear <- controls_vec_linear-1
    }
    names(controls_vec_linear) <- rownames(gep_mat)


    ratio_mat <- sapply(rev(sim_fracs), function(x){
      ratio_vec <- ((x*ctoi_vec_linear)+controls_vec_linear)/ctoi_vec_linear
      ratio_vec[is.infinite(ratio_vec)]
    })
    ratio_mat[is.nan(ratio_mat)] <- 1 # In case of 0/0


    ratio_vec <- controls_vec_linear/ctoi_vec_linear
    ratio_vec[is.nan(ratio_vec)] <- 1 # In case of 0/0
    ratio_vec[ratio_vec == Inf] <- 0

    ctoi_vec_linear["LYZ"]
    controls_vec_linear["LYZ"]
    controls_vec_linear["LYZ"]/ctoi_vec_linear["LYZ"]

    ctoi_vec_linear["MORC1"]
    controls_vec_linear["MORC1"]
    controls_vec_linear["MORC1"]/ctoi_vec_linear["MORC1"]

    ctoi_vec_linear["CXCR4"]
    controls_vec_linear["CXCR4"]
    controls_vec_linear["CXCR4"]/ctoi_vec_linear["CXCR4"]

    ctoi_vec_linear[names(ratio_vec[is.infinite(ratio_vec)])]
    controls_vec_linear["TLE6"]
    controls_vec_linear["GCK"]/ctoi_vec_linear["GCK"]

    ctoi_sim_list <- lapply(1:nSims, function(i){

      mix_vec <- Rfast::rowmeans(mix_linear[,sample(1:ncol(mix_linear), size = sample(2:ncol(mix_linear), 1), replace = FALSE)])
      mix_sub <- matrix(rep(mix_vec, length(sim_fracs)), ncol = length(sim_fracs), byrow = FALSE)
      rownames(mix_sub) <- rownames(mix)

      mix_sub["LYZ",] * ratio_mat["LYZ",]
      mix_sub["MORC1",] * ratio_mat["MORC1",]
      mix_sub["CXCR4",] * ratio_mat["CXCR4",]



      design_mat <- matrix(0, nrow = length(sim_fracs), ncol = length(sim_fracs))
      diag(design_mat) <- 1

      betamat <- matrix(rep(sizeEffect, length(sim_fracs)), ncol = length(sim_fracs), byrow = FALSE)
      rownames(betamat) <- rownames(mix)
      betamat <- log2((2^betamat-1) %*% diag(sim_fracs)+1)


      maxvec  <- apply(betamat, 1, max)


      simulation <- seqgendiff::thin_diff(mat = mix_sub, design_fixed = design_mat, coef_fixed = betamat)$mat
      rownames(simulation) <- rownames(mix)
      colnames(simulation) <- paste0("mix", "%%", sim_fracs)

      return(simulation)


      # good_genes <- unique(unlist(signatures[names(sort(c, decreasing = T)[1:5])]))
      # bad_genes <- unique(unlist(signatures[names(sort(c, decreasing = F)[1:5])]))
      # good_genes <- good_genes[!good_genes %in% bad_genes]
      # bad_genes <- bad_genes[!bad_genes %in% good_genes]
      #
      # cor(simulation[good_genes,ncol(simulation)], gep_mat[good_genes, ctoi])
      # cor(simulation[bad_genes,ncol(simulation)], gep_mat[bad_genes, ctoi])
      # sigs_genes <- unique(unlist(signatures[names(c)]))
      # cor(simulation[sigs_genes, ncol(simulation)], gep_mat[sigs_genes, ctoi])

      # sort(sapply(signatures_ctoi, function(sig){
      #   cor(simulation[sig,ncol(simulation)], ctoi_vec[sig])
      # }), decreasing = TRUE) -> xx
      # c[names(xx)]

      # all_cors <- apply(simulation, 2, function(mix){
      #   sapply(signatures_ctoi, function(sig){
      #     cor(mix[sig], ctoi_vec[sig])
      #   })
      # })

      # c[names(sort(rowMeans(all_cors), decreasing = TRUE))]


    })
    do.call(cbind, ctoi_sim_list)


  }, BPPARAM = param)
  names(sim_list) <- celltypes

  return(sim_list)




}
makeSimulations <- function(ref, labels, mix, gep_mat, cor_mat, dep_list, sim_fracs, sim_method, n_sims, ncores, seed2use){

  set.seed(seed2use)
  param <- BiocParallel::MulticoreParam(workers = ncores)


  celltypes <- unique(labels$label)

  getSubMatrix <- function(mat, sim_fracs, n_samples_sim){

    if (class(mat)[1] == "numeric") {
      ctoi_vec <- mat

    }else{
      ctoi_vec <- mat[,sample(1:ncol(mat), n_samples_sim)]
    }

    if (n_samples_sim > 1) {
      ctoi_vec <- Rfast::rowmeans(ctoi_vec)
      names(ctoi_vec) <- rownames(mat)
    }

    mat_sub <- matrix(rep(ctoi_vec, length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs))
    rownames(mat_sub) <- rownames(mat)
    return(mat_sub)
  }
  adjustLibSize <- function(ctoi_mat, controls_mat){

    # Scale to simulate counts data
    scale_factor <- 10000
    ctoi_mat_scaled <- round(ctoi_mat * scale_factor)
    controls_mat_scaled <- round(controls_mat * scale_factor)

    # Adjust reference-controls library size
    min_lib_size <- min(min(colSums(ctoi_mat_scaled)), min(colSums(controls_mat_scaled)))
    ref_ctoi_sub_lib_fracs <- min_lib_size/colSums(ctoi_mat_scaled)
    ref_controls_sub_lib_fracs <- min_lib_size/colSums(controls_mat_scaled)

    # Thin data to adjust lib size
    ctoi_mat_scaled_thin <- seqgendiff::thin_lib(ctoi_mat_scaled, thinlog2 = -log2(ref_ctoi_sub_lib_fracs), type = "thin")$mat
    controls_mat_scaled_thin <- seqgendiff::thin_lib(controls_mat_scaled, thinlog2 = -log2(ref_controls_sub_lib_fracs), type = "thin")$mat

    # Unscale
    ctoi_mat_thin <- ctoi_mat_scaled_thin/scale_factor
    rownames(ctoi_mat_thin) <- rownames(ctoi_mat)
    controls_mat_thin <- controls_mat_scaled_thin/scale_factor
    rownames(controls_mat_thin) <- rownames(controls_mat)

    return(list(ctoi_mat_thin = ctoi_mat_thin, controls_mat_thin = controls_mat_thin))
  }
  makeFractionMatrixControls <- function(controls_mat, sim_fracs, controls){


    # Multiple control mixture by fractions
    generateFractions <- function(target_sum, n_fracs = ncol(controls_mat)) {
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
        Rfast::rowmeans(controls_mat %*% diag(fracs))
      }
    })
    rownames(controls_frac_mat) <- rownames(ref)

    return(controls_frac_mat)
  }


  sim_list <- BiocParallel::bplapply(celltypes, function(ctoi){

    ref_ctoi <- ref[,labels$label == ctoi]

    # Get control cell types
    dep_cts <- unique(c(ctoi, unname(unlist(dep_list[[ctoi]]))))
    controls <- celltypes[!celltypes %in% dep_cts]


    # Generate n_sims simulations
    ctoi_sim_list <- lapply(1:n_sims, function(i){


      # Use n_samples_sim random samples for CTOI
      ref_ctoi_sub <- getSubMatrix(mat = ref_ctoi, sim_fracs, n_samples_sim)

        # Sample 1-10 control cell types
        controls2use <- sample(controls, sample(1:min(10, length(controls)), 1), replace = FALSE)
        # Generate controls matrix
        controls_mat <- sapply(controls2use, function(ctrl){

          control_samples <- labels[labels$label == ctrl,]$sample
          samples2use <- control_samples[1:sample(1:length(control_samples), 1)]

          if (length(samples2use) == 1) {
            ref[,samples2use]
          }else{
            Rfast::rowmeans(ref[,samples2use])
          }

        })
        rownames(controls_mat) <- rownames(ref)

        # Thin to adjust library size
        data_adjusted <- adjustLibSize(ctoi_mat = ref_ctoi_sub, controls_mat = controls_mat)
        ref_ctoi_sub <- data_adjusted$ctoi_mat_thin
        controls_mat <- data_adjusted$controls_mat_thin

        # Generete fractions ,at
        ctoi_frac_mat <- ref_ctoi_sub %*% diag(sim_fracs)
        control_frac_mat <- makeFractionMatrixControls(controls_mat, sim_fracs, controls = controls2use)


        simulation <- ctoi_frac_mat + control_frac_mat
        colnames(simulation) <- paste0("mix", "%%", sim_fracs)

        simulation



    })

    do.call(cbind, ctoi_sim_list)



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
      score <- singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
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

    xgb_params <- list(
      booster = "gbtree",  # Using linear booster
      alpha = 1,          # Lasso regularization term
      lambda = 1,            # No Ridge regularization
      eta = 0.01,             # Learning rate
      objective = "reg:squarederror",
      max_depth = 6,
      nthread = 1
    )
    model <- xgboost::xgboost(
      data = predictors,
      label = response,
      params = xgb_params,
      nrounds = 150,
      verbose = 0
    )
    #cor(round(predict(model, scores, type = "response"), 4), fracs, method = "spearman", use = "pairwise.complete.obs")

    importance_matrix <- xgboost::xgb.importance(feature_names = colnames(predictors), model = model)
    gains.tmp <- importance_matrix$Gain
    sigs_filtered <- c()
    for (i in 1:nrow(importance_matrix)) {
      p <- outliers::grubbs.test(gains.tmp)$p.value
      if (p >= 0.01) {
        break
      }
      sigs_filtered <- c(sigs_filtered, importance_matrix[i,1][[1]])
      gains.tmp <- gains.tmp[-1]
    }


    if (length(sigs_filtered) < 3) {
      sigs_filtered <- importance_matrix[1:3,1][[1]]
    }

    model <- xgboost::xgboost(
      data = predictors[,sigs_filtered],
      label = response,
      params = xgb_params,
      nrounds = 150,
      verbose = 0
    )
    # cor(round(predict(model, scores[,sigs_filtered], type = "response"), 4), fracs, method = "spearman", use = "pairwise.complete.obs")


    # Signature filtering with Lasso
    # cv_fit <- glmnet::cv.glmnet(data[,-ncol(data)], data[,ncol(data)], alpha = 1) # Using lasso penalty for feature selection
    # best_lambda <- cv_fit$lambda.min
    # coefficients <- as.matrix(coef(cv_fit, s = best_lambda))
    # sigs_filtered <- rownames(coefficients)[-1][coefficients[-1, , drop = FALSE] != 0]


    # RF <- RRF::RRF(x = as.data.frame(data[,-ncol(data)]), y = data[,ncol(data)], flagReg = 0, importance = TRUE, ntree = 1000)
    # RF_imp <- RF$importance[,"%IncMSE"] / max(RF$importance[,"%IncMSE"])
    # RRF <- RRF::RRF(x = data[,-ncol(data)], y = data[,ncol(data)], flagReg = 1, ntree = 1000, coefReg = (1-gamma) + gamma*RF_imp)
    # selected_features <- colnames(data)[RRF$feaSet]
    # data <- data[, c(selected_features, "frac")]

    # model <- randomForestSRC::var.select(frac ~ ., as.data.frame(data), method = "md", verbose = FALSE, refit = TRUE, fast = TRUE)
    # sigs_filtered <- model$topvars

    # model <- randomForestSRC::var.select(frac ~ ., as.data.frame(data)[,c(sigs_filtered, "frac")], method = "vh.vimp", verbose = FALSE, refit = TRUE, conservative = "high", fast = TRUE)
    # sigs_filtered <- model$topvars

    #cor(round(randomForestSRC::predict.rfsrc(model$rfsrc.refit.obj, newdata = as.data.frame(scores[,sigs_filtered]))$predicted, 4), fracs, method = "spearman", use = "pairwise.complete.obs")

    # Make sure there are at least three signatures
    # lasso_alpha <- 1
    # while(length(sigs_filtered) < 3) {
    #   lasso_alpha <- lasso_alpha - 0.1
    #   cv_fit <- glmnet::cv.glmnet(data[,-ncol(data)], data[,ncol(data)], alpha = lasso_alpha)
    #   coefficients <- as.matrix(coef(cv_fit, s = best_lambda))
    #   sigs_filtered <- rownames(coefficients)[-1][coefficients[-1, , drop = FALSE] != 0]
    # }

    # model <- glmnet::glmnet(data[,-ncol(data)], data[,ncol(data)], lambda = best_lambda, alpha = 0.5)
    # sigs_filtered <- colnames(data[,-ncol(data)])

    # Fit final model with Ridge
    # cv_fit <- glmnet::cv.glmnet(data[,sigs_filtered], data[,ncol(data)], alpha = 0)
    # best_lambda <- cv_fit$lambda.min
    # model <- glmnet::glmnet(data[,sigs_filtered], data[,ncol(data)], alpha = 0, lambda = best_lambda)

    # cor(round(predict(model, scores[,sigs_filtered], s = best_lambda, type = "response")[,1], 4), fracs, method = "spearman", use = "pairwise.complete.obs")

    # options(rf.cores=1, mc.cores=1)
    # model <- randomForestSRC::rfsrc.fast(frac ~ ., as.data.frame(data), forest = TRUE)
    # sigs_filtered <- colnames(data)[colnames(data) != "frac"]

    return(tibble(model = list(model), sigs_filtered = list(sigs_filtered)))

  }


  param <- BiocParallel::MulticoreParam(workers = ncores)

  #start <- Sys.time()
  models_list <- BiocParallel::bplapply(simulations_scored, function(data){
    fitModel(data)
  }, BPPARAM = param)
  #end <- Sys.time()
  #print(end-start)


  enframe(models_list, name = "celltype") %>%
    unnest(value) %>%
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
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, mix = NULL, ref_type, filtering_data = NULL, lineage_file = NULL, top_genes_frac = 1, medianGEP = TRUE, seed = 123, probs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4),
                        sim_fracs = c(0, seq(0.01, 0.25, 0.01), seq(0.3, 1, 0.05)), diff_vals = round(c(log2(1), log2(1.5), log2(2), log2(2.5), log2(3), log2(4), log2(5), log2(10), log2(20)), 3),
                        min_genes = 3, max_genes = 150, return_sigs = FALSE, return_sigs_filt = FALSE, sigsFile = NULL, minPBcells = 30, minPBsamples = 10,
                        ct_sims = 10, samples_frac = 0.1, simMethod = "ref_multi", nCores = 1){


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
    message("Filtering signatures by external datasets...")
    out <- filterSignatures(ref, labels, filtering_data, signatures, n_sigs = 3, ncores = nCores)
    # signatures <- out$filt_sigs
  }

  if (return_sigs_filt) {
    return(list(all_sigs = signatures,
         filt_sigs = out$filt_sigs,
         filt_cts = out$filt_cts))
  }

  # Make simulations
  message("Generating simulations...")
  simulations <- makeSimulations(labels, mix, gep_mat, sim_fracs, filtering_data, ncores = nCores, seed2use = seed)


  message("Scoring simulations...")
  # TODO: Fix problem when n_sims = 1
  simulations_scored <- scoreSimulations(signatures, simulations, n_sims = ncol(mix), nCores)


  # Filter signatures and train RF model
  message("Filtering signatures and training models...")
  models <- trainModels(simulations_scored, ncores = nCores, seed2use = seed)
  signatures_filt <- signatures[unlist(models$sigs_filtered)]
  models <- models[,-3]


  # Get spillover matrix
  # message("Generating spillover matrix...")
  # frac2use <- sim_fracs[which.min(abs(sim_fracs - 0.25))]
  # spill_mat <- getSpillOverMat(simulations, signatures_filt, dep_list, models, frac2use)


  # Save results in S4 object
  xCell2Sigs.S4 <- new("xCell2Signatures",
                       signatures = signatures_filt,
                       dependencies = dep_list,
                       models = models,
                       spill_mat = matrix(),
                       genes_used = rownames(ref)) ###################### Important!!!!


  message("Custom xCell2.0 reference ready!")

  return(xCell2Sigs.S4)

}
