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

# Generate a matrix of median expression of pure cell types
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

# This function return a correlation matrix given the counts and cell types
getCellTypeCorrelation <- function(pure_ct_mat, data_type){

  celltypes <- colnames(pure_ct_mat)

  # Use 50% most expressed genes because single-cell is zero inflated
  if (data_type != "sc") {
    mean_gene_expression <- Rfast::rowmeans(pure_ct_mat)
    high_gene_expression_cutoff <- quantile(mean_gene_expression, 0.5, na.rm=TRUE)
    top_expressed_gene <- mean_gene_expression > high_gene_expression_cutoff
    pure_ct_mat <- pure_ct_mat[top_expressed_gene,]
  }

  # Use top 85% most variable genes
  genes_var <- apply(pure_ct_mat, 1, var)
  most_var_genes_cutoff <- quantile(genes_var, 0.85, na.rm=TRUE)
  pure_ct_mat <- pure_ct_mat[genes_var > most_var_genes_cutoff,]

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

# This function return a vector of cell type dependencies
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


createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes){


  getSigs <- function(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes){

    # Remove dependent cell types
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, unname(unlist(dep_list[[type]])))]

    # Weight cell types by correlation
    type_weights <- cor_mat[type, not_dep_celltypes]

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

      # Score genes using weights
      gene_scores <- apply(diff_genes.mat, 1, function(x){
        sum(type_weights[which(x)])
      })

      gene_passed <- gene_scores[gene_scores > 0]

      # If less than min_genes passed move to next parameters
      if (length(gene_passed) < min_genes) {
        next
      }

      # Save signatures
      gene_passed <- sort(gene_passed, decreasing = TRUE)
      for (n_genes in round(seq(from = min_genes, to = max_genes, length.out = 8))) {

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
    getSigs(celltypes, type, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes)
  })
  all_sigs <- unlist(all_sigs)


  if (length(all_sigs) == 0) {
    warning("No signatures found for reference!")
  }

  # Make GeneSetCollection object
  signatures_collection <- GSEABase::GeneSetCollection(all_sigs)


  return(signatures_collection)
}


filterSignatures <- function(ref, labels, pure_ct_mat, dep_list, signatures_collection, mixture_fractions, grubbs_cutoff, simulations_cutoff){

  # This function created mixtures for the simulations
  getMixtures <- function(ref, labels, ct, ct_data, dep_list, max_control_type,
                          mixture_fractions){

    getControlsMeanExpression <- function(ref, labels, ct, dep_list, max_control_type, mixture_fractions){

      controls2use <- names(dep_list)[!names(dep_list) %in% c(ct, dep_list[[ct]])]

      n_ds <- length(unique(labels$dataset)) # This number of datasets will effect the sampling of controls
      n_slice <- ifelse(n_ds >= max_control_type, 1, round((max_control_type/n_ds)+0.5))
     if(n_ds == 1){
       n_slice <- length(controls2use)
      }
      if(length(controls2use) < max_control_type){
        max_control_type <- length(controls2use)
      }


      controls <- sapply(1:length(mixture_fractions), function(x){

        labels %>%
          filter(label %in% controls2use) %>% # Only controls
          group_by(dataset, label) %>%
          slice_sample(n=1) %>% # One cell type per dataset
          group_by(dataset) %>%
          slice_sample(n=n_slice) %>%
          ungroup() %>%
          slice_sample(n=max_control_type) %>%
          pull(sample) %>%
          ref[,.] %>%
          as.matrix() %>%
          Rfast::rowmeans()


      })

      controls_fracs <- controls %*% diag(1-mixture_fractions)
      return(controls_fracs)

    }

    mixSmaples <- function(ref, labels, ct, dep_list, max_control_type, ct_mean_expression, mixture_fractions){

      # Generate a matrix of CTOI fractions:
      m <- matrix(rep(ct_mean_expression, length(mixture_fractions)), byrow = FALSE, ncol = length(mixture_fractions)) %*% diag(mixture_fractions)
      # Get random controls fractions
      c <- getControlsMeanExpression(ref, labels, ct, dep_list, max_control_type, mixture_fractions)
      # Add controls
      m <- m + c

      rownames(m) <- rownames(ref)
      colnames(m) <- as.character(mixture_fractions)

      return(m)

    }


    ct_data %>%
      group_by(dataset) %>%
      summarise(samples = list(sample)) %>%
      rowwise() %>%
      mutate(n_samples = length(samples)) %>%
      ungroup() %>%
      # Use top 20 references with most samples
      top_n(n = 20, wt = n_samples) %>%
      rowwise() %>%
      # Get cell type mean expression for each dataset
      mutate(ct_mean_expression = list(Rfast::rowmeans(as.matrix(ref[,samples])))) %>%
      mutate(mixtures = list(mixSmaples(ref, labels, ct, dep_list, max_control_type, ct_mean_expression, mixture_fractions))) %>%
      pull(mixtures) %>%
      return()
  }


  getMixturesCors <- function(signatures_collection, ct, mixture_ranked, mixture_fractions){

    sigs <- signatures_collection[startsWith(names(signatures_collection), paste0(ct, "#"))]

    # Score every ranked mixtures of CTOI
    cors <- sapply(mixture_ranked, function(ranked_mix){

      # Score
      scores <- sapply(sigs, simplify = TRUE, function(sig){
        singscore::simpleScore(ranked_mix, upSet = sig, centerScore = FALSE)$TotalScore
      })
      if (is.list(scores)) {
        sigs <- sigs[-which(lengths(scores) == 0)]
        scores <- sapply(sigs, simplify = TRUE, function(sig){
          singscore::simpleScore(ranked_mix, upSet = sig, centerScore = FALSE)$TotalScore
        })
      }
      colnames(scores) <- names(sigs)

      # Correlation
      apply(scores, 2, function(sig_scores){
        cor(sig_scores, mixture_fractions, method = "spearman")
      })

    })

    median_cors <- Rfast::rowMedians(cors)
    names(median_cors) <- rownames(cors)

    return(median_cors)
  }




  scores_mat <- matrix(nrow = length(signatures_collection),
                       ncol = ncol(pure_ct_mat),
                       dimnames = list(names(signatures_collection), colnames(pure_ct_mat)))

  sig_type <- unlist(lapply(strsplit(names(signatures_collection), "#"), "[", 1))

  for (type in unique(sig_type)) {
    print(type)

    type_signatures <- signatures_collection[type == sig_type]


    dep_cells <- unname(unlist(dep_list[[type]]))

    if (length(dep_cells) > 0) {
      types_to_use <- !colnames(scores_mat) %in% dep_cells
    }else{
      types_to_use <- rep(TRUE, ncol(scores_mat))
    }

    sub_mix <- pure_ct_mat[,types_to_use]
    sub_mix_ranked <- singscore::rankGenes(sub_mix)
    for (i in 1:length(type_signatures)) {
      sig <- type_signatures[i]
      scores_out <- singscore::simpleScore(sub_mix_ranked, upSet = sig[[1]], centerScore = FALSE)$TotalScore
      scores_mat[which(rownames(scores_mat) == names(sig)), colnames(sub_mix_ranked)] <- scores_out
    }

  }

  scores_mat_tidy <- scores_mat %>%
    as_tibble(., rownames = NA) %>%
    rownames_to_column(var = "signature") %>%
    pivot_longer(cols = -signature, values_to = "score", names_to = "sample_ct") %>%
    separate(signature, into = "signature_ct", sep = "#", remove = FALSE, extra = "drop")%>%
    drop_na()

  # Filter by simulations
  simulations <- labels %>%
    group_by(ont, label) %>%
    nest() %>%
    rowwise() %>%
    # TODO: takes too long
    mutate(mixture = list(getMixtures(ref = ref, labels = labels, ct = label, ct_data = data, dep_list = dep_list, max_control_type = 5, mixture_fractions = mixture_fractions))) %>%
    mutate(mixture_ranked = list(lapply(mixture, function(mix){singscore::rankGenes(mix)})))

  simulations.cors <- simulations %>%
    mutate(mixture_cor = list(getMixturesCors(signatures_collection = signatures_collection, ct = label, mixture_ranked = mixture_ranked, mixture_fractions = mixture_fractions)))

  simulations.filtered <- simulations.cors %>%
    mutate(sig_filtered = list(names(mixture_cor)[mixture_cor >= quantile(mixture_cor, simulations_cutoff)])) %>%
    pull(sig_filtered) %>%
    unlist()


  # Filter by Grubb's test
  grubbs <- scores_mat_tidy %>%
   group_by(signature_ct, signature) %>%
   summarise(grubbs_statistic = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$statistic[1]) %>%
   filter(grubbs_statistic >= quantile(grubbs_statistic, grubbs_cutoff)) %>%
   pull(signature)



  signatures_collection_filtered <- signatures_collection[names(signatures_collection) %in% simulations.filtered]

  filter_signature_out <- list("scoreMatTidy" = scores_mat_tidy, "sigCollectionFilt" = signatures_collection_filtered)

  return(filter_signature_out)
}


#' @slot labels ...
#' @slot dependencies ...
#' @slot all_signatures ...
#' @slot filtered_signatures ...
#' @name xCell2Signatures
#' @importFrom methods new
# Create S4 object for the new reference
setClass("xCell2Signatures", slots = list(
  labels = "data.frame",
  dependencies = "list",
  all_signatures = "GeneSetCollection",
  filtered_signatures = "GeneSetCollection"
))

# Remove - for debugging
if (0 == 1) {
  data_type = "sc"; lineage_file = NULL; mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02))
  probs = c(.1, .25, .33333333, .5); diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5); min_genes = 5; max_genes = 200
}


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
#' @importFrom xCell2 xCell2GetLineage
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
#' @return An S4 object containing the signatures, cell type labels, and cell type dependencies.
#' @export
xCell2Train <- function(ref, labels, data_type, lineage_file = NULL, mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02)),
                        probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                        min_genes = 5, max_genes = 500){


  # Validate inputs
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
    dep_list <- xCell2::xCell2GetLineage(labels = labels[,1:2], out_file = NULL)
  }else{
    dep_list <- getDependencies(lineage_file)
  }

  # Generate signatures for each cell type
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs, dep_list, include_descendants = TRUE)
  message("Generating signatures...")
  signatures_collection <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes)

  # Filter signatures
  message("Filtering signatures...")
  filter_signature_out <- filterSignatures(ref, labels, pure_ct_mat, dep_list, signatures_collection, mixture_fractions, grubbs_cutoff = 0.8, simulations_cutoff = 0.8)
  scores_mat_pure_tidy <- filter_signature_out$scoreMatTidy
  signatures_collection_filtered <- filter_signature_out$sigCollectionFilt


  # TODO: Weight signatures with Elastic Net


  # TODO: Linear transformation


  # TODO: Spillover correction


  xCell2Ref.S4 <- new("xCell2Signatures", labels = labels, dependencies = dep_list, all_signatures = signatures_collection, filtered_signatures = signatures_collection_filtered)

  return(xCell2Ref.S4)

}
