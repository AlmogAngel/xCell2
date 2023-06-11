# Signatures generation analysis for scRNA-seq benchmarking

library(tidyverse)
library(xCell2)
data("ts_labels_with_ontology")

# Load dataset
datasets <- sub(".rds", "", list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/"))
ds <- "Thymus"

ds.rds <- readRDS(paste0("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/", ds, ".rds"))
ref <- ds.rds$singleCellExpr
labels <- tibble(ont = "NA",
                 label = ds.rds$singleCellLabels,
                 sample = colnames(ref),
                 dataset = ds.rds$singleCellSubjects) %>%
  mutate(ont = plyr::mapvalues(sample, ts_labels_with_ontology$sample, ts_labels_with_ontology$ont.fine, warn_missing = FALSE),
         label =  gsub("_", "-", label))
labels <- as.data.frame(labels)
bulk <- ds.rds$bulk
truth <- ds.rds$bulkRatio




# Generate signatures for dataset:
# (!) Run first functions of xCell2 for: dep_list, quantiles_matrix and cor_mat

# Method A
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
signatures_collection_A <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes = TRUE)

#Method B
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
      diff <- param.df[i, ]$diff_vals # difference criteria
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
signatures_collection_B <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes, weight_genes = FALSE)



# Get signatures correlations
getSigsCors <- function(sigs, bulk, truth, labels, method_info){

  bulk_ranked <- singscore::rankGenes(bulk)
  ds_celltypes <- unique(labels$label)


  # Score signatures
  ds_scores_list <- pbapply::pbsapply(ds_celltypes, function(ctoi){
    signatures_ctoi <- sigs[startsWith(names(sigs), paste0(ctoi, "#"))]
    scores <- t(sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    }))
    colnames(scores) <- colnames(bulk_ranked)
    rownames(scores) <- names(signatures_ctoi)
    scores
  })

  # Get correlations
  rownames(truth) <- gsub("_", "-", rownames(truth))
  ds_cors_list <- pbapply::pbsapply(names(ds_scores_list), function(ctoi){
    ds_ctoi_sigs_scores <- ds_scores_list[[ctoi]]

    sigs_cors <- apply(ds_ctoi_sigs_scores, 1, function(sig_scores){
      cor(sig_scores, truth[ctoi,], method = "spearman")
    })
  })

  # Make tidy
  enframe(ds_cors_list, name = "celltype") %>%
    unnest_longer(value, indices_to = "signature", values_to = "cor") %>%
    mutate(method = method_info) %>%
    return(.)
}


cor_A <- getSigsCors(sigs = signatures_collection_A, bulk, truth, labels, method_info = "weights")

cor_B <- getSigsCors(sigs = signatures_collection_B, bulk, truth, labels, method_info = "no_weights")


# Plot correlations
rbind(cor_A, cor_B) %>%
  ggplot(., aes(x=celltype, y=cor, fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "green", size=0.5)



cor_A %>%
  filter(celltype == ctoi) %>%
  separate(signature, into = c("name", "prob", "diff", "n_genes"), sep = "_") %>%
  ggplot(., aes(x=diff, y=cor, fill=diff)) +
  geom_boxplot()
