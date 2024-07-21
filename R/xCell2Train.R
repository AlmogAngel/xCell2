validateInputs <- function(ref, labels, ref_type){

  if (length(unique(labels$label)) < 3) {
    stop("Reference must have at least 3 cell types!")
  }

  if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
    stop("Reference must be one of those classes: matrix, dgCMatrix, Matrix")
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
      tmp <- ref[,ctoi_samples]
      colnames(tmp) <- as.character(1:ncol(tmp))
      tmp
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
prepRefMix <- function(ref, mix, ref_type, min_sc_genes, human2mouse){

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

    if (min(ref) < 3) {
      # Adding 3 to reference restrict inclusion of small changes
      ref <- ref + 3
    }

  }else{
    if(max(ref) < 50){

      ref <- 2^ref
      if (min(ref) == 1) {
        ref <- ref-1
      }

      if (min(ref) < 3) {
        # Adding 3 to reference restrict inclusion of small changes
        ref <- ref + 3
      }

    }else{
      if (min(ref) < 3) {
        # Adding 3 to reference restrict inclusion of small changes
        ref <- ref + 3
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
createSignatures <- function(labels, dep_list, quantiles_matrix, probs, diff_vals, min_genes, max_genes, min_frac_ct_passed, num_threads){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, diff_vals, min_genes, max_genes, min_frac_ct_passed){

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

    type.sigs <- getSigs(celltypes, type, dep_list, quantiles_matrix, probs, diff_vals,
                         min_genes, max_genes, min_frac_ct_passed)


    # Check for minimum 3 signatures per cell type
    if (length(type.sigs) < 3) {
      # Relax prob-diff-min_genes parameters
      probs2use <- c(probs, max(probs)*1.25, max(probs)*1.5, max(probs)*2)
      probs2use <- probs2use[probs2use < 0.5]
      diff_vals2use <- c(min(diff_vals)*0, min(diff_vals)*0.5, min(diff_vals)*0.75, diff_vals)
      min_genes2use <- round(min_genes*0.5)
      min_genes2use <- ifelse(min_genes2use < 3, 3, min_genes2use)

      type.sigs <- getSigs(celltypes, type, dep_list, quantiles_matrix, probs = probs2use, diff_vals = diff_vals2use,
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
learnParams <- function(gep_mat, cor_mat, signatures, dep_list, ref_type, top_spill_value, sc_spill_relaxing_factor, num_threads){

  param <- BiocParallel::MulticoreParam(workers = num_threads)

  celltypes <- colnames(gep_mat)
  gep_mat_linear <- 2^gep_mat
  sim_fracs <- c(0, seq(0.01, 0.25, 0.01))
  frac2use <- 0.25

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
      controls <- celltypes[celltypes != ctoi]
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
      controls <- celltypes[celltypes != ctoi]
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
    if (!is.null(dep_list)) {
      dep_cts <- dep_cts[dep_cts != ctoi]
      final_scores[dep_cts] <- 0
    }
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
#' @importFrom Rfast rowMedians rowmeans rowsums Sort
#' @importFrom Matrix rowMeans rowSums colSums Diagonal
#' @importFrom singscore rankGenes simpleScore
#' @param ref A reference gene expression matrix (genes in rows samples/cells in columns).
#' @param mix description
#' @param labels A data frame in which the rows correspond to samples in the ref. The data frame must have four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545" or NA if there is no ontology).
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample that match the column name in ref.
#'   "dataset": sample's source dataset or subject (for single-cell).
#' @param ref_type Gene expression data type: "rnaseq" for bulk RNA-Seq, "array" for micro-array, or "sc" for scRNA-Seq.
#' @param seed Set seed for reproducible results (optional).
#' @param min_pb_cells For scRNA-Seq reference only - minimum number of cells in the pseudo-bulk (optional).
#' @param min_pb_samples For scRNA-Seq reference only - minimum number of pseudo-bulk samples (optional).
#' @param min_sc_genes Minimum number of genes for scRNA-Seq (default 10000).
#' @param use_ontology A Boolean for using ontological integration (TRUE)
#' @param lineage_file Path to the cell type lineage file generated with `xCell2GetLineage` function (optional).
#' @param num_threads Number of threads for parallel processing.
#' @param human2mouse A Boolean for converting human genes to mouse genes.
#' @param top_spill_value Maximum spillover compensation correction value
#' @param sc_spill_relaxing_factor description
#' @param return_signatures A Boolean to return just the signatures.
#' @param return_analysis A Boolean to return the xCell2Analysis results (do not return signatures object).
#' @param use_sillover A Boolean to use spillover correction in xCell2Analysis (return_analysis much be TRUE)
#' @param spillover_alpha A numeric for spillover alpha value in xCell2Analysis.
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref,
                        mix = NULL,
                        labels,
                        ref_type,
                        human2mouse = FALSE,
                        lineage_file = NULL,
                        seed = 123,
                        num_threads = 1,
                        use_ontology = TRUE,
                        return_signatures = FALSE,
                        return_analysis = FALSE,
                        use_sillover = TRUE,
                        spillover_alpha = 0.25,
                        min_pb_cells = 30,
                        min_pb_samples = 10,
                        min_sc_genes = 1e4,
                        top_spill_value = 0.5,
                        sc_spill_relaxing_factor = 0.1
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

  # Prepare reference and mixture data: human to mouse genes transformation, normalization,log2 transformation, shared genes
  out <- prepRefMix(ref, mix, ref_type, min_sc_genes, human2mouse)
  ref <- out$ref.out
  mix <- out$mix.out
  shared_genes <- rownames(ref)

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
  probs <- c(0.1, 0.25, 0.333, 0.49)
  diff_vals <- c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5)
  min_genes <- 8
  max_genes <- 200
  min_frac_ct_passed <- 0.5
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs, num_threads)
  message("Generating signatures...")
  signatures <- createSignatures(labels, dep_list, quantiles_matrix, probs, diff_vals, min_genes, max_genes, min_frac_ct_passed, num_threads)

  if (return_signatures) {
    xCell2.S4 <- new("xCell2Object",
                     signatures = signatures,
                     dependencies = list(),
                     params = data.frame(),
                     spill_mat = matrix(),
                     genes_used = shared_genes)
    return(xCell2.S4)
  }

  # Learn linear transformation parameters
  message("Learning linear transformation and spillover parameters...")
  gep_mat <- makeGEPMat(ref, labels)
  cor_mat <- getCellTypeCorrelation(gep_mat, ref_type)
  params <- learnParams(gep_mat, cor_mat, signatures, dep_list, ref_type, top_spill_value, sc_spill_relaxing_factor, num_threads)


  # Save results in S4 object
  if (is.null(dep_list)) {
    dep_list <- list()
  }
  xCell2.S4 <- new("xCell2Object",
                   signatures = signatures,
                   dependencies = dep_list,
                   params = params$params,
                   spill_mat = params$spillmat,
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
