###

library(tidyverse)

getTopVariableGenes <- function(ref, min_genes, sensitivity){

  ref.srt <- Seurat::CreateSeuratObject(counts = ref)
  ref.srt <- Seurat::FindVariableFeatures(ref.srt, selection.method = "vst")
  plot1 <- Seurat::VariableFeaturePlot(ref.srt)
  genesVar <- plot1$data$variance.standardized
  names(genesVar) <- rownames(plot1$data)
  genesVar <- sort(genesVar, decreasing = TRUE)
  idx <- 1:length(genesVar)
  KneeCutOff <- kneedle::kneedle(idx, genesVar, sensitivity = sensitivity)[2]
  #ggplot(data.frame(x=idx, y=genesVar), aes(x=x, y=y)) +
  #  geom_point()+
  #  geom_hline(yintercept = KneeCutOff, col="red")
  topGenesVar <- names(genesVar[genesVar >= KneeCutOff])

  if (length(topGenesVar) < min_genes) {
    topGenesVar <- names(genesVar[1:min_genes])
  }
  return(topGenesVar)
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

  if (data_type != "sc") {
    # Use top 3K highly variable genes
    mean_gene_expression <- Rfast::rowmeans(pure_ct_mat)
    high_gene_expression_cutoff <- quantile(mean_gene_expression, 0.5, na.rm=TRUE) # Cutoff for top 50% expression genes
    top_expressed_gene <- mean_gene_expression > high_gene_expression_cutoff
    genes_sd <- apply(pure_ct_mat[top_expressed_gene,], 1, sd)
    sd_cutoff <-  sort(genes_sd, decreasing = TRUE)[1001]# Get top 1K genes with high SD as highly variable genes
    pure_ct_mat <- pure_ct_mat[genes_sd > sd_cutoff,]
  }

  # Make correlation matrix
  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

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
    dep_cells <- c(strsplit(pull(ont[i,3]), ";")[[1]], strsplit(pull(ont[i,4]), ";")[[1]])
    dep_cells <- gsub("_", "-", dep_cells)
    dep_list[[i]] <- dep_cells[!is.na(dep_cells)]
  }

  return(dep_list)
}

# Generate a list with quantiles matrices for each cell type
makeQuantiles <- function(ref, labels, probs){

  celltypes <- unique(labels[,2])
  samples <- labels[,2]

  quantiles_matrix <- lapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    # If there is one sample for this cell type - duplicate the sample to make a data frame
    if (sum(type_samples) == 1) {
      type.df <- cbind(ref[,type_samples], ref[,type_samples])
    }else{
      type.df <- ref[,type_samples]
    }
    quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1-probs))), na.rm=TRUE))
  })
  names(quantiles_matrix) <- celltypes

  return(quantiles_matrix)
}



createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes){

  celltypes <- unique(labels[,2])
  samples <- labels[,2]

  all_sigs <- list()
  for (type in celltypes){
    # Get dependent cell types and remove them from the quantiles matrix.
    dep_cells <- dep_list[[type]]
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, dep_cells)]

    # Find similar cell types
    cor_values <- cor_mat[type, ]
    cor_values <- sort(cor_values[!names(cor_values) %in% c(type, dep_cells)], decreasing = TRUE)

    # Remove cell-types with dependencies in cor_values (new!)
    to_remove <- c()
    for (i in 1:length(cor_values)) {
      if (names(cor_values[i]) %in% to_remove) {
        next
      }
      deps <- dep_list[[names(cor_values[i])]]
      to_remove <- c(to_remove, names(cor_values)[names(cor_values) %in% deps])
    }
    cor_values <- cor_values[!names(cor_values) %in% to_remove]

    # TODO: Should we add extra genes for similar cell-types or not? What are similar cell-types?
    cor_cells_cutoff <- quantile(cor_values, 0.9, na.rm=TRUE)
    sim_cells <- names(cor_values[which(cor_values > cor_cells_cutoff & cor_values > 0.85)])
    # sim_cells <- names(cor_values[which(cor_values > 0.5)])
    # sim_cells <- c()

    for (diff in diff_vals) {
      for (p in 1:length(probs)) {

        # Get a list of Boolean matrices with genes that pass the quantiles criteria
        diff_genes <- lapply(not_dep_celltypes, function(x){
          get(type, quantiles_matrix)[p,] > get(x, quantiles_matrix)[nrow(quantiles_matrix[[1]])-p+1,] + diff
        })

        diff_genes.mat <- matrix(unlist(diff_genes), nrow = length(diff_genes), byrow = TRUE,
                                 dimnames = list(not_dep_celltypes, rownames(ref)))


        # Find signature genes for similar cell types
        sim_genes <- c()
        n_sim <- 0

        # In case there is only one similar cell type
        if (length(sim_cells) == 1) {
          n_sim <- 1
          sim_genes <- names(which(diff_genes.mat[rownames(diff_genes.mat) == sim_cells,] > 0))
          if (length(sim_genes) > round(max_genes*0.25)) {
            sim_genes <- sim_genes[1:round(max_genes*0.25)]
          }
        }

        # In case there are more than one cell types
        if (length(sim_cells) > 1) {
          for (n_sim in length(sim_cells):1) {
            sim_genes <- names(which(colSums(diff_genes.mat[rownames(diff_genes.mat) %in% sim_cells,]) >= n_sim))
            if (length(sim_genes) > 0) {
              # Make sure sim_genes is not bigger than 1/4 of max_genes
              if (length(sim_genes) > round(max_genes*0.25)) {
                sim_genes <- sim_genes[1:round(max_genes*0.25)]
              }
              break
            }
          }
        }


        # Find signature genes for all cell types
        for (n_all in nrow(diff_genes.mat):round(nrow(diff_genes.mat)*0.5)) {
          genes <- names(which(colSums(diff_genes.mat) >= n_all))
          # Merge with sim_genes
          genes <- unique(c(sim_genes, genes))
          # If check there are enough genes
          if (length(genes) < min_genes) {
            next
          }
          # Check if there are too many genes
          if (length(genes) > max_genes) {
            genes <- genes[1:max_genes]
          }
          # Save signature
          if (length(genes) > 0) {
            sig_name <-  paste(paste0(type, "#"), probs[p], diff, n_sim, n_all, sep = "_")
            all_sigs[[sig_name]] <- GSEABase::GeneSet(genes, setName = sig_name)
          }
        }

      }

    }

  }

  if (length(all_sigs) == 0) {
    warning("No signatures found for reference!")
  }

  # Make GeneSetCollection object
  signatures_collection <- GSEABase::GeneSetCollection(all_sigs)

  # Check if every cell type got at least one signature
  sig_type <- unlist(lapply(strsplit(names(signatures_collection), "#"), "[", 1))
  no_sig_cts <- celltypes[!celltypes %in% sig_type]

  if (length(no_sig_cts) != 0) {
    warning(paste("No signatures found for cell type(s): ", no_sig_cts, collapse = " "))
  }

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

        # labels %>%
        #   filter(label %in% controls2use) %>%
        #   group_by(dataset) %>%
        #   slice_sample(n=1) %>% # One cell type per dataset
        #   group_by(label) %>%
        #   slice_sample(n=1) %>% # One sample per datasets
        #   ungroup() %>%
        #   slice_sample(n=max_control_type) %>%
        #   pull(sample) %>%
        #   ref[,.] %>%
        #   as.matrix() %>%
        #   Rfast::rowmeans()

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

    print(ct) #remove

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



  celltypes <- colnames(pure_ct_mat)

  scores_mat <- matrix(nrow = length(signatures_collection),
                       ncol = ncol(pure_ct_mat),
                       dimnames = list(names(signatures_collection), colnames(pure_ct_mat)))

  sig_type <- unlist(lapply(strsplit(names(signatures_collection), "#"), "[", 1))

  for (type in celltypes) {

    if (sum(type == sig_type) == 0) {
      errorCondition(paste0("No signatures found for cell type: ", type))
    }

    type_signatures <- signatures_collection[type == sig_type]
    dep_cells <- dep_list[[type]]

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
    mutate(mixture = list(getMixtures(ref = ref, labels = labels, ct = label, ct_data = data, dep_list = dep_list, max_control_type = 5, mixture_fractions = mixture_fractions))) %>%
    mutate(mixture_ranked = list(lapply(mixture, function(mix){singscore::rankGenes(mix)})))

  simulations.cors <- simulations %>%
    mutate(mixture_cor = list(getMixturesCors(signatures_collection = signatures_collection, ct = label, mixture_ranked = mixture_ranked, mixture_fractions = mixture_fractions)))

  simulations.filtered <- simulations.cors %>%
    mutate(sig_filtered = list(names(mixture_cor)[mixture_cor >= quantile(mixture_cor, simulations_cutoff)])) %>%
    pull(sig_filtered) %>%
    unlist()


  # Filter by Grubb's test
  #grubbs <- scores_mat_tidy %>%
  #  group_by(signature_ct, signature) %>%
  #  summarise(grubbs_statistic = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$statistic[1]) %>%
  #  filter(grubbs_statistic >= quantile(grubbs_statistic, grubbs_cutoff)) %>%
  #  pull(signature)



  signatures_collection_filtered <- signatures_collection[names(signatures_collection) %in% simulations.filtered]

  filter_signature_out <- list("scoreMatTidy" = scores_mat_tidy, "sigCollectionFilt" = signatures_collection_filtered)

  return(filter_signature_out)
}

