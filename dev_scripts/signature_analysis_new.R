library(xCell2)
library(tidyverse)
library(ggpubr)
library(Seurat)

setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")

# Function to load validation dataset
loadVals <- function(valList,
                     valMixDir = "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/",
                     valTruthDir = "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/",
                     conversion_mat = celltype_conversion){

  # load mixtures
  valMixtures <- lapply(valList, function(tissue){
    mixtures <- lapply(tissue, function(val){
      as.matrix(read.table(paste0(valMixDir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
    })
    names(mixtures) <- tissue
    mixtures
  })

  # Load truth
  valTruth <- lapply(valList, function(tissue){
    truths <- lapply(tissue, function(val){
      truth <- as.matrix(read.table(paste0(valTruthDir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
      # Change cell types labels
      rownames(truth) <- plyr::mapvalues(rownames(truth), conversion_mat$all_labels, conversion_mat$xCell2_labels, warn_missing = FALSE)
      rows <- rownames(truth)
      truth <- suppressWarnings(apply(truth, 2, as.numeric))
      rownames(truth) <- rows
      truth[is.na(truth)] <- 0
      truth <- truth[rowSums(truth) != 0,]
      rows <- rownames(truth)
      truth[!duplicated(rows),]
    })
    names(truths) <- tissue
    truths
  })

  return(list(mixtures = valMixtures, truth = valTruth))
}

# Function to combine reference-validation data
combineRefVal <- function(vals, refList, splitDependencies){


  vals.tbl <- enframe(vals$mixtures, name = "val_type", value = "mixture") %>%
    unnest_longer(mixture, indices_to = "val_dataset") %>%
    left_join(unnest_longer(enframe(vals$truth, name = "val_type", value = "truth"), truth, indices_to = "val_dataset"), by = c("val_type", "val_dataset"))

  refs.tbl <- enframe(refList, name = "ref_type", value = "ref_name") %>%
    unnest_longer(ref_name, indices_to = "ref_tissue")

  types <- unique(vals.tbl$val_type)

  # Generate a tibble with validations-references combinations
  combined.tbl <- lapply(types, function(type){
    val_tmp <- vals.tbl[vals.tbl$val_type %in% type,]
    type <- if(type == "other") "mixed" else c(type, "mixed")
    ref_tmp <- refs.tbl[refs.tbl$ref_tissue %in% type,]
    crossing(val_tmp, ref_tmp)
  }) %>%
    do.call(rbind, .) %>%
    dplyr::select(-mixture, -truth)
  rm(vals.tbl, refs.tbl)


  # Load refs dependencies
  refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
  getDependencies <- function(lineage_file_checked){
    ont <- readr::read_tsv(lineage_file_checked, show_col_types = FALSE) %>%
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
  refsDepList <- lapply(refList, function(ref_type){
    refsDep <- lapply(ref_type, function(ref){

      # Load reference
      ref <- readRDS(paste0("references/", ref, "_ref.rds"))

      # Get cell types dependencies
      dep_list <- getDependencies(ref$lineage_file)

      dep_cts <- sapply(names(dep_list), function(d){
        deps <- unname(unlist(dep_list[[d]]))
      })
      dep_cts[order(sapply(dep_cts, length), decreasing = T)]
    })
    names(refsDep) <- ref_type
    refsDep
  })


  # Split references with dependent cell types for deconvolution methods
  if (!splitDependencies) {

    combined.tbl %>%
      rowwise() %>%
      # Get number of samples in the validation dataset
      mutate(n_val_samples = ncol(vals$truth[[val_type]][[val_dataset]])) %>%
      # Get shared cell types between the validation and reference
      mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), names(refsDepList[[ref_type]][[ref_name]])))) %>%
      mutate(n_shared_celltypes = length(shared_celltypes)) %>%
      # Must have at least 3 cell types
      filter(n_shared_celltypes > 2) %>%
      mutate(refsDeps = NA) %>%
      # I non deconvolution method celltype_classes = shared_celltypes
      mutate(celltype_classes = list(shared_celltypes)) %>%
      return()

  }else{

    splitDepCellTypes <- function(refsDeps){

      refsDeps <- refsDeps[order(sapply(refsDeps, length), decreasing = T)]

      classes <- list("1" = c())
      for (ct in names(refsDeps)) {

        classes <- classes[order(sapply(classes, length), decreasing = F)]
        deps <- refsDeps[[ct]]

        for (class in 1:length(classes)) {
          if (!any(deps %in% classes[[class]])) {
            classes[[class]] <- c(classes[[class]], ct)
            break
          }
        }

        if (!ct %in% unname(unlist(classes))) {
          classes[[paste0(length(classes)+1)]] <- ct
        }
      }

      one_element_list <- unname(lapply(classes, length) == 1)
      if (sum(one_element_list) > 0) {
        classes[!one_element_list][[1]] <- unname(c(unlist(classes[!one_element_list][1]), unlist(classes[one_element_list])))
        classes <- classes[!one_element_list]
      }


      return(classes)
    }

    combined.tbl %>%
      rowwise() %>%
      # Get number of samples in the validation dataset
      mutate(n_val_samples = ncol(vals$truth[[val_type]][[val_dataset]])) %>%
      # Get shared cell types between the validation and reference
      mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), names(refsDepList[[ref_type]][[ref_name]])))) %>%
      mutate(n_shared_celltypes = length(shared_celltypes)) %>%
      # Must have at least 3 cell types
      filter(n_shared_celltypes > 2) %>%
      # Create a list of dependent cell types - filter cell types that are shared between the reference and validation
      mutate(refsDeps = list(refsDepList[[ref_type]][[ref_name]][shared_celltypes])) %>%
      # Split dependent cell types in each CIBERSORTx run
      mutate(celltype_classes = list(splitDepCellTypes(refsDeps))) %>%
      unnest(celltype_classes) %>%
      return()

  }

}


# Function to get correlations from all unfiltered signatures
getCors <- function(vals, refsRDSList, shared_celltypes, valType, valName, refType, refName){

  val_ref <- paste0(valName, "_", refName)
  print(val_ref)

  mix <- vals$mixtures[[valType]][[valName]]
  ref <- refsRDSList[[refType]][[refName]]$ref
  labels <- refsRDSList[[refType]][[refName]]$labels
  lineage_file <- refsRDSList[[refType]][[refName]]$lineage_file

  if (refType == "sc") {
    shared_cleaned_genes <- xCell2CleanGenes(ref = ref, mix = mix, top_var_genes = TRUE, use_protein_coding = TRUE, n_var_genes = 10000)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <- xCell2CleanGenes(ref = ref, mix = mix, use_protein_coding = FALSE)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }

  sigs.s4 <- readRDS(paste0("/bigdata/almogangel/xCell2_data/dev_data/", val_ref, "_sigs.rds"))
  sigs <- sigs.s4@signatures

  if (!all(sigs.s4@genes_used == rownames(ref))) {
    stop()
  }

  truth <- cyto.vals$truth[[valType]][[valName]]

  shared_samples <- intersect(colnames(mix), colnames(truth))
  truth <- truth[,shared_samples]
  mix <- mix[,shared_samples]

  signatures_collection_tmp <- sigs[gsub("#.*", "", names(sigs)) %in% rownames(truth)]

  mix_ranked <- singscore::rankGenes(mix)

  sigs.cor.mat <- sapply(1:length(signatures_collection_tmp), function(i){

    sig <- signatures_collection_tmp[i]
    celltype <- gsub("#.*", "",  names(sig))
    scores <- suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig[[1]], centerScore = FALSE))
    samples <- rownames(scores)
    scores <- scores[,1]
    names(scores) <- samples
    truth_ct <- truth[celltype,]

    if (!all(names(scores) == names(truth_ct))) {
      print("Error!! mix != truth - samples names")
    }

    if (length(scores) == 0) {
      cor.out <- NA
    }else{
      cor.out <- cor(scores, truth_ct, method = "spearman")
    }
    c(names(sig), celltype, cor.out)
  }) %>%
    t()
  colnames(sigs.cor.mat) <- c("signature", "celltype", "cor")

  as_tibble(sigs.cor.mat) %>%
    mutate(cor = as.numeric(cor)) %>%
    return(.)

}


# Load cell type conversion
celltype_conversion <- read_tsv("celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# Pick references to use
refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood",
                       tumor = "sc_pan_cancer"))

# Pick validation datasets
cyto.vals.list <- list(blood = c("BG_blood", "GSE107011", "GSE107572", "GSE127813", "GSE53655", "GSE60424"),
                       tumor = c("ccRCC_cytof_CD45+", "NSCLC_cytof", "WU_ccRCC_RCCTC"),
                       other = c("GSE120444", "GSE115823"))

# Load validation datasets
cyto.vals <- loadVals(cyto.vals.list)

# Combine reference-validation data
refval.tbl <- combineRefVal(vals = cyto.vals, refList, splitDependencies = FALSE)




# (1) Generate signatures form every reference ------------

# print("Generating signatures...")
#
# refs <- c("bp", "kass_tumor", "kass_blood", "lm22", "sc_pan_cancer", "ts_blood")
#
# for (ref in refs) {
#   print(ref)
#
#   if (ref %in% c("bp", "kass_tumor", "kass_blood")) {
#     dataType <- "rnaseq"
#   }else if(ref == "lm22"){
#     dataType <- "array"
#   }else if(ref %in% c("sc_pan_cancer", "ts_blood")){
#     dataType <- "sc"
#   }
#
#   ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
#   all_sigs <- xCell2Train(ref.in$ref, ref.in$labels, data_type = dataType, lineage_file = ref.in$lineage_file, filter_sigs = FALSE)
#   saveRDS(all_sigs, paste0("/bigdata/almogangel/xCell2_data/dev_data/", ref, "_all_sigs.rds"))
# }
# print("Done - all sigs")



# (2) Calculate correlations ------------


print("Calculating correlation...")

refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
refsRDSList <- lapply(refList, function(ref_type){
  refs <- lapply(ref_type, function(ref){
    # Load reference
    ref.in <- readRDS(paste0("references/", ref, "_ref.rds"))
    ref.in
  })
  names(refs) <- ref_type
  refs
})

# i = 1
# valType=x[i,]$val_type[[1]]; valName=x[i,]$val_dataset[[1]]; refType=x[i,]$ref_type[[1]]; refName=x[i,]$ref_name[[1]]; shared_celltypes=x[i,]$shared_celltypes[[1]]


# Get all signatures correlations with reference-validation data combination
all.cyto.cors <- refval.tbl %>%
  dplyr::select(val_type:shared_celltypes) %>%
  rowwise() %>%
  mutate(cors_data = list(getCors(cyto.vals, refsRDSList, shared_celltypes, valType = val_type, valName = val_dataset, refType = ref_type, refName = ref_name)))

all.cyto.cors %>%
  ungroup() %>%
  unnest(cors_data) %>%
  saveRDS(., "/bigdata/almogangel/xCell2_data/dev_data/all.cyto.cors.rds")




# (3) Analyze diff/prob/n_genes ------------

cors.tbl <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/all.cyto.cors.rds")

# Diff vs. Prob correlation heatmap
cors.tbl %>%
  filter(ref_type == "sc") %>%
  separate(signature, into = c("removeMe", "prob", "diff", "n_genes"), sep = "_", remove = FALSE) %>%
  dplyr::select(-n_genes, -signature, -removeMe) %>%
  group_by(ref_name, val_dataset, prob, diff) %>%
  summarise(cor = mean(cor)) %>%
  group_by(ref_name, prob, diff) %>%
  summarise(cor = mean(cor)) %>%
  group_by(prob, diff) %>%
  summarise(cor = mean(cor, na.rm = T)) %>%
  ggplot(., aes(x = prob, y = diff, fill = cor)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(x = "Prob", y = "Diff", fill = "Cor")

cors.tbl %>%
  filter(ref_type != "sc") %>%
  separate(signature, into = c("removeMe", "prob", "diff", "n_genes"), sep = "_", remove = FALSE) %>%
  dplyr::select(-n_genes, -signature, -removeMe) %>%
  group_by(ref_name, val_dataset, prob, diff) %>%
  summarise(cor = mean(cor)) %>%
  group_by(ref_name, prob, diff) %>%
  summarise(cor = mean(cor)) %>%
  group_by(prob, diff) %>%
  summarise(cor = mean(cor, na.rm = T)) %>%
  ggplot(., aes(x = prob, y = diff, fill = cor)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(x = "Prob", y = "Diff", fill = "Cor")


cors.tbl %>%
  separate(signature, into = c("removeMe", "prob", "diff", "n_genes"), sep = "_", remove = FALSE) %>%
  dplyr::select(-prob, -diff, -signature, -removeMe) %>%
  group_by(ref_name, val_dataset, n_genes) %>%
  summarise(cor = mean(cor)) %>%
  group_by(ref_name, n_genes) %>%
  summarise(cor = mean(cor)) %>%
  group_by(n_genes) %>%
  summarise(cor = mean(cor, na.rm = T)) %>%
  mutate(n_genes = as.numeric(n_genes)) %>%
  mutate(n_genes_cat = cut(n_genes,
                           breaks = seq(0, max(n_genes), by = 10),
                           include.lowest = TRUE,
                           labels = paste(seq(0, max(n_genes)-10, by = 10),
                                          seq(10, max(n_genes), by = 10),
                                          sep = "-"))) %>%
  ggplot(., aes(x = n_genes_cat, y = cor, fill=n_genes_cat)) +
  geom_boxplot()



# (4) Analyze filtering ------------


# Generate and score simulations

# Run xCell2 first steps
refs <- c("bp", "kass_tumor", "kass_blood", "lm22", "sc_pan_cancer", "ts_blood")
refs.data <- list()
for (ref_name in refs) {
  print(ref_name)
  data_type <- ifelse(ref_name %in% c("sc_pan_cancer", "ts_blood"), "sc",
                      ifelse(ref_name == "lm22", "array", "rnaseq"))

  ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref_name, "_ref.rds"))
  signatures_collection <- readRDS(paste0("/bigdata/almogangel/xCell2_data/dev_data/", ref_name, "_all_sigs.rds"))


  ref = ref.in$ref
  labels = ref.in$labels
  lineage_file = ref.in$lineage_file
  clean_genes = TRUE
  diff_vals = c(1, 1.32, 1.585, 2, 3, 4, 5)
  min_genes = 5
  max_genes = 200


  # Load functions from xCell2Train
  validateInputs <- function(ref, labels, data_type){
    if (length(unique(labels$label)) < 3) {
      stop("Reference must have at least 3 cell types")
    }

    if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
      stop("ref must be one of those classes: matrix, dgCMatrix, Matrix")
    }

    if (!"data.frame" %in% class(labels)) {
      stop("labels must be a dataframe.")
    }

    if (!data_type %in% c("rnaseq", "array", "sc")) {
      stop("data_type should be 'rnaseq', 'array' or 'sc'.")
    }

    if (sum(grepl("_", labels$label)) != 0) {
      message("Changing underscores to dashes in cell-types labels!")
      labels$label <- gsub("_", "-", labels$label)
    }

    out <- list(ref = ref,
                labels = labels)
    return(out)

  }
  cleanGenes <- function(ref, gene_groups){
    gene.list <- hs.genelist
    if (all(startsWith(rownames(ref), "ENSG"))) {
      message("Cleaning genes from reference (assuming genes are in Ensembl ID)...")
      ref <- ref[!rownames(ref) %in% gene.list[gene.list$gene_group %in% gene_groups, 2],]
      return(ref)
    }else{
      message("Cleaning genes from reference (assuming genes are in symbol ID)...")
      ref <- ref[!rownames(ref) %in% gene.list[gene.list$gene_group %in% gene_groups, 3],]
      return(ref)
    }
  }
  NormalizeRef <- function(ref, data_type){

    if (data_type == "sc") {

      # Log-Normalize
      message("Normalizing and transforming scRNA-Seq reference to log1p-space.")

      # TODO: Change log1p to log2(x+1)
      genes_names <- rownames(ref)
      ref.srt <- CreateSeuratObject(counts = ref)
      ref.srt <- NormalizeData(ref.srt, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
      ref.norm <-  ref.srt@assays$RNA@data
      rownames(ref.norm) <- genes_names # Because Seurat change genes names from "_" to "-"

      return(ref.norm)

    }else{

      if(max(ref) >= 50){
        message("Transforming reference to log2-space (maximum expression value >= 50).")
        ref.norm <- log2(ref+1)
        return(ref.norm)
      }else{
        message("Assuming reference is already in log2-space (maximum expression value < 50).")
        return(ref)
      }

    }

  }
  makePureCTMat <- function(ref, labels, use_median){

    celltypes <- unique(labels$label)

    pure_ct_mat <- sapply(celltypes, function(type){
      type_samples <- labels[,2] == type
      if (sum(type_samples) == 1) {
        type_vec <- as.vector(ref[,type_samples])
      }else{
        if(use_median){
          type_vec <- if("matrix" %in% class(ref)) Rfast::rowMedians(ref[,type_samples]) else sparseMatrixStats::rowMedians(ref[,type_samples])
        }else{
          type_vec <- if("matrix" %in% class(ref)) Rfast::rowmeans(ref[,type_samples]) else Matrix::rowMeans(ref[,type_samples])
        }
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


  # Validate inputs
  inputs_validated <- validateInputs(ref, labels, data_type)
  ref <- inputs_validated$ref
  labels <- inputs_validated$labels

  # Clean genes
  if (clean_genes) {
    gene_groups <- c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY")
    ref <- cleanGenes(ref, gene_groups)
  }

  # Normalize Reference
  ref <- NormalizeRef(ref, data_type)

  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  pure_ct_mat <- makePureCTMat(ref, labels, use_median = TRUE)
  cor_mat <- getCellTypeCorrelation(pure_ct_mat, data_type)

  # Get cell type dependencies list
  message("Loading dependencies...")
  if (is.null(lineage_file)) {
    dep_list <- xCell2::xCell2GetLineage(labels, out_file = NULL)
  }else{
    dep_list <- getDependencies(lineage_file)
  }


  # Store ref data
  refs.data[[ref_name]] <- list("pure_ct_mat" = pure_ct_mat, "cor_mat" = cor_mat,
                                "dep_list" = dep_list, "signatures_collection" = signatures_collection)
}

sim_fracs = c(0, 0.001, 0.005, seq(0.01, 0.25, 0.03))
makeSimulations <- function(pure_ct_mat, cor_mat, dep_list, sim_fracs, simple_sim){

  celltypes <- colnames(pure_ct_mat)

  if (simple_sim) {
    # Simple simulations: (1) use highest fraction. (2) use pure cell type expression matrix. (3) use one control.
    sim_frac <- sim_fracs[length(sim_fracs)]

    # Make CTOIs fractions matrix
    ctoi_mat <- pure_ct_mat * sim_frac

    # Find a control for each cell type
    controls <- unname(sapply(celltypes, function(ctoi){
      dep_cts <- unname(unlist(dep_list[[ctoi]]))
      not_dep_cts <- celltypes[!celltypes %in% dep_cts]
      names(sort(cor_mat[ctoi, not_dep_cts])[1])
    }))

    # Make controls fractions matrix
    controls_mat <- sapply(controls, function(ctrl){
      pure_ct_mat[,ctrl] * (1-sim_frac)
    })

    # Combine CTOI and control matrices
    sim_names <- paste0(colnames(ctoi_mat), "%%", colnames(controls_mat))
    sim_mat <- ctoi_mat + controls_mat
    colnames(sim_mat) <- sim_names

    # # In case there is one control for all other cell types -> make a second controls matrix just for him
    # controls_abundance <- sort(table(controls), decreasing = TRUE)
    # if(controls_abundance[1] == length(celltypes)-1){
    #   abundant_control <- names(controls_abundance[1])
    #   controls2 <- unname(sapply(celltypes, function(ctoi){
    #     dep_cts <- unname(unlist(dep_list[[ctoi]]))
    #     not_dep_cts <- celltypes[!celltypes %in% dep_cts]
    #     not_dep_cts <- not_dep_cts[not_dep_cts != abundant_control]
    #     names(sort(cor_mat[ctoi, not_dep_cts])[1])
    #   }))
    #   controls_mat2 <- sapply(controls2, function(ctrl){
    #     pure_ct_mat[,ctrl] * (1-sim_frac)
    #   })
    #
    #   sim_names2 <- paste0(colnames(ctoi_mat), "%%", colnames(controls_mat2))
    #   sim_mat2 <- ctoi_mat + controls_mat2
    #   colnames(sim_mat2) <- sim_names2
    #
    #   sim_list <- list(sim1 = list(sim_mat = sim_mat, controls_mat = controls_mat), sim_mat2 = sim_mat2)
    #
    # }else{
    #   sim_list <- list(sim1 = list(sim_mat = sim_mat, controls_mat = controls_mat), sim_mat2 = NULL)
    # }


    return(sim_mat)
  }



}
filterSignatures <- function(sim = simple_sim_mat, signatures_collection, dep_list, n_null = 1000, n_cpu = 10){

  scoreSim <- function(sim, signatures_collection, dep_list, n_null){

    # Rank mixture
    sim_ranked <- singscore::rankGenes(sim)

    # Generate null distribution for p-values
    sigs_genes <- unique(unlist(signatures_collection))
    non_sigs_genes <- rownames(sim_ranked)[!rownames(sim_ranked) %in% sigs_genes]
    sigs_lengths <- unique(lengths(signatures_collection))
    all_lengths_null_scores <- pbapply::pblapply(sigs_lengths, function(len){
      tmp_genes <- sample(non_sigs_genes, len)
      singscore::generateNull(
        upSet = tmp_genes,
        rankData = sim_ranked,
        subSamples = 1:ncol(sim_ranked),
        centerScore = FALSE,
        B = n_null,
        ncores = n_cpu,
        seed = 1)
    })
    names(all_lengths_null_scores) <- sigs_lengths


    # Score signatures and get p-values
    sigs_scores_pvals <- pbapply::pblapply(signatures_collection, function(sig){
      sig_len <- length(sig)
      scoredf <- singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)
      pvals <- singscore::getPvals(all_lengths_null_scores[[as.character(sig_len)]], scoredf)
      tibble("sim_name" = names(pvals), "score" = scoredf$TotalScore, "pval" = pvals)
    })
    names(sigs_scores_pvals) <- names(signatures_collection)


    # Make results tidy
    sigs_scores_pvals_tidy <- enframe(sigs_scores_pvals, name = "signature") %>%
      unnest(value) %>%
      separate(sim_name, into = c("sim_celltype", "sim_control"), sep = "%%") %>%
      separate(signature, into = "sig_celltype", sep = "#", extra = "drop", remove = FALSE)

    # Clean scores
    sigs_scores_pvals_tidy_clean <- sigs_scores_pvals_tidy %>%
      filter(sig_celltype != sim_control) %>% # Signature cell type cannot be the same as the control
      rowwise() %>%
      filter(!sim_celltype %in% unname(unlist(dep_list[[sig_celltype]])) & !sim_control %in% unname(unlist(dep_list[[sig_celltype]])))  # sim_celltype and sim_control cannot be dependent on signature cell type

    return(sigs_scores_pvals_tidy_clean)

  }


  sigs_scores_pvals_tidy_clean <- scoreSim(sim, signatures_collection, dep_list, n_null)

  # Get top 25% or top 10 signatures by p-value delta
  top_pval_sigs <- sigs_scores_pvals_tidy_clean %>%
    mutate(pval = -log(pval)) %>%
    mutate(pval_type = ifelse(sig_celltype == sim_celltype, "ctoi", "controls")) %>%
    group_by(sig_celltype, signature, pval_type) %>%
    summarise(pval = mean(pval)) %>%
    pivot_wider(names_from = pval_type, values_from = pval) %>%
    mutate(delta_pval = `ctoi` - `controls`) %>%
    group_by(sig_celltype) %>%
    top_n(n = max(10, 0.25*n()), wt = delta_pval) %>%
    pull(signature)


  # Filter by top_pval_sigs and get top 10% or top 5 signatures by p-value delta
  filtered_sigs <- sigs_scores_pvals_tidy_clean %>%
    filter(signature %in% top_pval_sigs) %>%
    mutate(score_type = ifelse(sig_celltype == sim_celltype, "ctoi", "controls")) %>%
    group_by(sig_celltype, signature, score_type) %>%
    summarise(score = mean(score)) %>%
    pivot_wider(names_from = score_type, values_from = score) %>%
    mutate(delta_score = `ctoi` - `controls`) %>%
    group_by(sig_celltype) %>%
    top_n(n = max(5, 0.1*n()), wt = delta_score) %>%
    pull(signature)


  # celltypes <- unique(gsub("#.*", "", names(signatures_collection)))
  # table(gsub("#.*", "", filtered_sigs))
  # table(gsub("#.*", "", unique(sigs_scores_pvals_tidy_clean$signature)))

  signatures_collection_filtered <- signatures_collection[names(signatures_collection) %in% filtered_sigs]

  out <- list(sim_results = sigs_scores_pvals_tidy_clean,
              filtered_sigs = signatures_collection_filtered)

  return(out)
}

all.scores.pvals.list <- pbapply::pblapply(refs, function(ref_name){
  print(ref_name)
  simple_sim_mat2 <- makeSimulations(pure_ct_mat = refs.data[[ref_name]]$pure_ct_mat, cor_mat = refs.data[[ref_name]]$cor_mat,
                                     dep_list = refs.data[[ref_name]]$dep_list, sim_fracs, simple_sim = TRUE)
  out <- filterSignatures(sim = simple_sim_mat2, signatures_collection = refs.data[[ref_name]]$signatures_collection,
                          dep_list = refs.data[[ref_name]]$dep_list)
  out$sim_results <- out$sim_results %>%
    mutate(ref = ref_name)
  out
})
names(all.scores.pvals.list) <- refs
saveRDS(all.scores.pvals.list, "/bigdata/almogangel/xCell2_data/dev_data/all.scores.pvals.list.rds")



# Check current filtering criteria

# Load correlations
ref.cors.tbl <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/all.cyto.cors.rds") %>%
  mutate(ref_sig = paste0(ref_name, "#", signature),
         filter = "before")


top_ref_sig <- ref.cors.tbl %>%
  group_by(ref_name, celltype, signature) %>%
  summarise(cor = cor) %>%
  group_by(ref_name, celltype) %>%
  top_n(n = max(10, 0.1*n()), wt = cor) %>%
  mutate(ref_sig = paste0(ref_name, "#", signature)) %>%
  pull(ref_sig)


all.scores.pvals.list <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/all.scores.pvals.list.rds")
tmp <- lapply(all.scores.pvals.list, function(x){names(x$filtered_sigs)})
ref_sig_filtered <- unlist(lapply(names(tmp), function(x){paste0(x, "#", tmp[[x]])}))


ref.cors.tbl %>%
  filter(ref_sig %in% ref_sig_filtered) %>%
  mutate(filter = "after") %>%
  rbind(., ref.cors.tbl) %>%
  #filter(ref_name == "bp") %>%
  ggplot(., aes(x=filter, y=cor, fill=filter)) +
  geom_boxplot()



# Try new filtering criteria


all_sim_results <- lapply(all.scores.pvals.list, function(x){x$sim_results}) %>%
  do.call(rbind, .) %>%
  ungroup() %>%
  mutate(sigs = "all sigs")

all_sim_results %>%
  mutate(ref_sig = paste0(ref, "#", signature)) %>%
  filter(ref_sig %in% top_ref_sig) %>%
  dplyr::select(-ref_sig) %>%
  mutate(sigs = "top sig") %>%
  rbind(all_sim_results, .) %>%
  filter(sig_celltype == sim_celltype) %>%
  #mutate(pval = -log(pval)) %>%
  ggplot(., aes(x=ref, y=score, fill=sigs)) +
  geom_boxplot()


y <- all_sim_results %>%
  mutate(pval = -log(pval)) %>%
  mutate(pval_type = ifelse(sig_celltype == sim_celltype, "ctoi", "controls")) %>%
  group_by(ref, sig_celltype, signature, pval_type) %>%
  summarise(pval = mean(pval)) %>%
  pivot_wider(names_from = pval_type, values_from = pval) %>%
  mutate(delta_pval = `ctoi` - `controls`) %>%
  mutate(ref_sig = paste0(ref, "#", signature)) %>%
  mutate(sigs = "all sig")


y %>%
  filter(ref_sig %in% top_ref_sig) %>%
  mutate(sigs = "top sig") %>%
  rbind(y, .) %>%
  #mutate(pval = -log(pval)) %>%
  ggplot(., aes(x=ref, y=delta_pval, fill=sigs)) +
  geom_boxplot()


# Get top signatures by p-value delta
top_pval_sigs <- all_sim_results %>%
  mutate(pval = -log(pval)) %>%
  mutate(pval_type = ifelse(sig_celltype == sim_celltype, "ctoi", "controls")) %>% View()
group_by(ref, sig_celltype, signature, pval_type) %>%
  summarise(pval = mean(pval)) %>%
  pivot_wider(names_from = pval_type, values_from = pval) %>%
  mutate(delta_pval = `ctoi` - `controls`) %>%
  group_by(ref, sig_celltype) %>%
  top_n(n = max(10, 0.25*n()), wt = delta_pval) %>%
  #top_frac(0.25, wt = delta_pval) %>%
  mutate(ref_sig = paste0(ref, "#", signature)) %>%
  pull(ref_sig)


# Filter by top_pval_sigs and get top  signatures by scores delta
ref_sig_filtered <- all_sim_results %>%
  mutate(ref_sig = paste0(ref, "#", signature)) %>%
  #filter(ref_sig %in% top_pval_sigs) %>%
  mutate(score_type = ifelse(sig_celltype == sim_celltype, "ctoi", "controls")) %>%
  group_by(ref, sig_celltype, signature, score_type) %>%
  summarise(score = mean(score)) %>%
  pivot_wider(names_from = score_type, values_from = score) %>%
  mutate(delta_score = `ctoi` - `controls`) %>%
  group_by(ref, sig_celltype) %>%
  top_n(n = max(5, 0.1*n()), wt = delta_score) %>%
  #top_frac(0.1, wt = delta_score) %>%
  mutate(ref_sig = paste0(ref, "#", signature)) %>%
  pull(ref_sig)


ref.mean.cors.tbl %>%
  filter(ref_sig %in% ref_sig_filtered) %>%
  mutate(filter = "after") %>%
  rbind(., ref.mean.cors.tbl) %>%
  group_by(filter) %>%
  summarise(cor = median(cor, na.rm = T))

ggplot(., aes(x=filter, y=cor, fill=filter)) +
  geom_boxplot()










after.filtering <- before.filtering %>%
  filter(ref_sig %in% filtered_sigs) %>%
  mutate(filter = "Top 25% p-values + Top 10% delta_score")


after.filtering %>%
  rbind(., before.filtering) %>%
  ggplot(., aes(x=ref, y=cor, fill=filter)) +
  geom_boxplot()


























# old -----------------------




# Delta score
delta_score <- before.filtering %>%
  mutate(scores_type = ifelse(sig_celltype == mix_celltype, "ctoi", "controls")) %>%
  group_by(ref, signature, scores_type, cor) %>%
  summarise(score = mean(score)) %>%
  pivot_wider(names_from = scores_type, values_from = score) %>%
  mutate(delta_score = `ctoi` - `controls`) %>%
  mutate(filter = "before")


delta_score %>%
  ggplot(., aes(x = delta_score, y = cor)) +
  geom_point(aes(col=ref)) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  stat_cor(method = "pearson", cor.coef.name = "r", label.x = 0.1, label.y = -0.6) +
  labs(x = "delta_score", y = "cor") +
  theme_minimal()


delta_score %>%
  filter(signature %in% pvalues_sigs) %>%
  group_by(ref) %>%
  top_frac(0.1, wt=delta_score) %>%
  mutate(filter = "Top 25% p-values + Top 10% delta_score") %>%
  rbind(., delta_score) %>%
  #filter(ref %in% c("kass_blood", "ts_blood")) %>%
  ggplot(., aes(x=filter, y=cor, fill=filter)) +
  geom_boxplot()

delta_score_sigs <- delta_score %>%
  group_by(ref) %>%
  top_frac(0.1, wt=delta_score) %>%
  mutate(filter = "Top 10% delta_score") %>%
  pull(signature)


# p-value
pvalues_max_filter_sigs <- before.filtering %>%
  mutate(pval = -log(pval)) %>%
  mutate(pval_type = ifelse(sig_celltype == mix_celltype, "ctoi", "controls")) %>%
  group_by(ref, signature, pval_type, cor) %>%
  summarise(pval = max(pval)) %>%
  pivot_wider(names_from = pval_type, values_from = pval) %>%
  mutate(delta_pval = `ctoi` - `controls`) %>%
  filter(delta_pval > 0) %>%
  pull(signature)

pvalues <- before.filtering %>%
  filter(signature %in% pvalues_max_filter_sigs) %>%
  mutate(pval_type = ifelse(sig_celltype == mix_celltype, "ctoi", "controls")) %>%
  group_by(ref, signature, pval_type, cor) %>%
  summarise(pval = mean(pval)) %>%
  mutate(pval = -log(pval)) %>%
  pivot_wider(names_from = pval_type, values_from = pval) %>%
  mutate(delta_pval= `ctoi` - `controls`) %>%
  mutate(filter = "before")


pvalues %>%
  filter(delta_pval > 0) %>%
  ggplot(., aes(x = delta_pval, y = cor)) +
  geom_point(aes(col=ref)) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  stat_cor(method = "pearson", cor.coef.name = "r", label.x = 0.1, label.y = -0.6) +
  labs(x = "delta_pval", y = "cor") +
  theme_minimal()


pvalues %>%
  group_by(ref) %>%
  top_frac(0.1, wt=delta_pval) %>%
  mutate(filter = "Top 10% delta_pval") %>%
  rbind(., pvalues) %>%
  ggplot(., aes(x=filter, y=cor, fill=filter)) +
  geom_boxplot()


pvalues %>%
  #filter(signature %in% delta_score_sigs) %>%
  group_by(ref) %>%
  top_frac(0.1, wt=delta_pval) %>%
  mutate(filter = "Top 10% delta_pval") %>%
  rbind(., pvalues) %>%
  ggplot(., aes(x=filter, y=cor, fill=filter)) +
  geom_boxplot()

pvalues_sigs <- pvalues %>%
  group_by(ref) %>%
  top_frac(0.25, wt=delta_pval) %>%
  mutate(filter = "Top 25% delta_pval") %>%
  pull(signature)




mixture_fractions <- c(0, 0.001, 0.005, seq(0.01, 0.25, 0.03))
mix_frac <- mixture_fractions[length(mixture_fractions)] # Use the highest fraction


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

  # Rank mixture
  mix_ranked <- singscore::rankGenes(mix)

  # Generate null distribution for p-values
  sigs_genes <- unique(unlist(signatures_collection))
  non_sigs_genes <- rownames(mix_ranked)[!rownames(mix_ranked) %in% sigs_genes]
  sigs_lengths <- unique(lengths(signatures_collection))
  all_lengths_null_scores <- pbapply::pblapply(sigs_lengths, function(len){
    tmp_genes <- sample(non_sigs_genes, len)
    singscore::generateNull(
      upSet = tmp_genes,
      rankData = mix_ranked,
      subSamples = 1:ncol(mix_ranked),
      centerScore = FALSE,
      B = 1000,
      ncores = 40,
      seed = 1)
  })
  names(all_lengths_null_scores) <- sigs_lengths


  sigs_scores_pvals <- pbapply::pblapply(signatures_collection, function(sig){
    sig_len <- length(sig)
    scoredf <- singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)
    pvals <- singscore::getPvals(all_lengths_null_scores[[as.character(sig_len)]], scoredf)
    tibble("mix_name" = names(pvals), "score" = scoredf$TotalScore, "pval" = pvals)
  })
  names(sigs_scores_pvals) <- names(signatures_collection)

  sigs_scores_pvals_tidy <- enframe(sigs_scores_pvals, name = "signature") %>%
    unnest(value) %>%
    separate(mix_name, into = c("mix_celltype", "mix_control"), sep = "%%") %>%
    separate(signature, into = "sig_celltype", sep = "#", extra = "drop", remove = FALSE)


  # Clean scores
  sigs_scores_pvals_tidy_clean <- sigs_scores_pvals_tidy %>%
    filter(sig_celltype != mix_control) %>% # Signature cell type cannot be the same as the control
    rowwise() %>%
    filter(!mix_celltype %in% unname(unlist(dep_list[[sig_celltype]])) & !mix_control %in% unname(unlist(dep_list[[sig_celltype]])))  # Simulation CTOI/control cannot be dependent on signature cell type

  return(sigs_scores_pvals_tidy_clean)

}


all.scores.pvals.list <- pbapply::pblapply(refs, function(ref){
  print(ref)
  mix_list <- makeMixture(pure_ct_mat = refs.data[[ref]]$pure_ct_mat, cor_mat = refs.data[[ref]]$cor_mat,
                          dep_list = refs.data[[ref]]$dep_list, mix_frac)
  sigs_scores_pvals_tidy_clean <- scoreMixture(mix = mix_list$mix1$mix_mat, signatures_collection = refs.data[[ref]]$signatures_collection,
                                               dep_list = refs.data[[ref]]$dep_list)
  sigs_scores_pvals_tidy_clean %>%
    mutate(ref = ref)
})
names(all.scores.pvals.list) <- refs

saveRDS(all.scores.pvals.list, "/bigdata/almogangel/xCell2_data/dev_data/all.cors.pvals.list.rds")



#

all.cors.pvals.list <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/all.cors.pvals.list.rds")
all.cors.pvals.tbl <- all.cors.pvals.list %>%
  do.call(rbind, .)


cors.tbl <- cors.tbl %>%
  mutate(ref_sig = paste0(ref, "#", signature))

all.cors.pvals.tbl <- all.cors.pvals.tbl %>%
  ungroup() %>%
  mutate(ref_sig = paste0(ref, "#", signature))


all.cors.pvals.tbl <- all.cors.pvals.tbl %>%
  filter(ref_sig %in% cors.tbl$ref_sig)


length(unique(all.cors.pvals.tbl$ref_sig))
length(unique(cors.tbl$ref_sig))
length(intersect(unique(cors.tbl$ref_sig), unique(all.cors.pvals.tbl$ref_sig)))

ref_sigs_cors <- cors.tbl %>%
  group_by(ref_sig) %>%
  summarise(cor = mean(cor))




# (1) First filtering - Top of delta score between CTOI and median score (of all other cell types)

all.cors.pvals.tbl <- all.cors.pvals.tbl %>%
  ungroup() %>%
  mutate(logpVal = -log(pval))



# mean_sig_score <- all.cors.pvals.tbl %>%
#   ungroup() %>%
#   filter(sig_celltype != mix_celltype) %>%
#   group_by(ref_sig) %>%
#   summarise(mean_score = mean(score))
#
# mean_sig_pval <- all.cors.pvals.tbl %>%
#   ungroup() %>%
#   filter(sig_celltype != mix_celltype) %>%
#   group_by(ref_sig) %>%
#   summarise(mean_pval= mean(logpVal))
#
# all.scores.pvals <- all.cors.pvals.tbl %>%
#   ungroup() %>%
#   filter(sig_celltype == mix_celltype) %>%
#   rowwise() %>%
#   left_join(mean_sig_score, by = "ref_sig") %>%
#   left_join(mean_sig_pval, by = "ref_sig") %>%
#   mutate(delta_score = score - mean_score) %>%
#   mutate(delta_pval = logpVal - mean_pval) %>%
#   ungroup()
#
#
# final <- all.scores.pvals %>%
#   left_join(ref_sigs_cors, by="ref_sig")

final <- all.cors.pvals.tbl %>%
  ungroup() %>%
  filter(sig_celltype == mix_celltype) %>%
  left_join(ref_sigs_cors, by="ref_sig")


# cut_points <- seq(min(final$logpVal, na.rm = TRUE) %/% 2 * 2, max(final$logpVal, na.rm = TRUE) %/% 2 * 2 + 2, by = 2)
# cut_labels <- paste(cut_points[-length(cut_points)], cut_points[-1]-1, sep = "-")
#
# # Cut logpVal into groups based on the defined cut points and labels
# final$logpVal_group <- cut(final$logpVal, breaks = cut_points, include.lowest = TRUE, labels = cut_labels)

# Create the plot
final %>%
  group_by(ref) %>%
  mutate(is_top10 = cor >= quantile(cor, 0.90, na.rm = TRUE)) %>%
  ggplot(., aes(x = ref, y = logpVal, fill = is_top10)) +
  geom_boxplot()


final %>%
  ggplot(., aes(x=delta_pval, y=cor)) +
  geom_point()

ggplot(final, aes(x = logpVal, y = cor)) +
  geom_point(aes(col=ref)) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  stat_cor(method = "pearson", cor.coef.name = "r", label.x = 0.1, label.y = -0.6) +
  labs(x = "logpVal", y = "cor") +
  theme_minimal()



write_csv(final, file = "/bigdata/almogangel/xCell2_data/dev_data/final.tsv")

sigs_filt1 <- all.cors.pvals.tbl %>%
  group_by(ref, sig_celltype) %>%
  top_frac(n=0.1, wt=delta_score) %>%
  pull(signature)


sigs_filt2 <- all.cors.pvals.tbl %>%
  group_by(ref, sig_celltype) %>%
  top_frac(n=0.1, wt=delta_pval) %>%
  pull(signature)


filt1 <- cors.tbl %>%
  filter(signature %in% sigs_filt1) %>%
  mutate(filt = "delta_score")

filt2 <- cors.tbl %>%
  filter(signature %in% sigs_filt2) %>%
  mutate(filt = "delta_pval")

filt3 <- cors.tbl %>%
  filter(ref == "bp") %>%
  filter(signature %in% sigs_filt3) %>%
  mutate(filt = "simulations")

cors.tbl %>%
  filter(ref == "bp") %>%
  mutate(filt = "before") %>%
  rbind(., filt3) %>%
  mutate(filt = factor(filt, level=c("before", "simulations"))) %>%
  ggplot(., aes(x=filt, y=cor)) +
  geom_boxplot(aes(fill=filt), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 45, hjust=1),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),) +
  labs(x = "", y = "Spearman r", x = NULL,  fill = "filter")


# Simulations .............

mixture_fractions <- c(0, 0.001, 0.005, seq(0.01, 0.25, 0.03))
ref = "bp"
signatures_collection <- readRDS(paste0("/bigdata/almogangel/xCell2_data/dev_data/", ref, "_all_sigs.rds"))
ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
ref <- ref.in$ref
labels <- ref.in$labels
dep_list <- getDependencies(ref.in$lineage_file)
pure_ct_mat <- makePureCTMat(ref.in$ref, ref.in$labels, use_median = TRUE)
cor_mat <- getCellTypeCorrelation(pure_ct_mat, data_type)


sim_list <- makeSimulations(ref, labels, mixture_fractions, dep_list, cor_mat, n_ct_sim = 10, add_noise = FALSE)
scores_list <- scoreCTOISimulations(signatures = signatures_collection, sim_list)

scores_all_sims_tidy <- enframe(scores_list, name = "celltype") %>%
  unnest_longer(value, indices_to = "sim_id", values_to = "scores") %>%
  rowwise() %>%
  mutate(scores = list(mat2tidy(scores))) %>%
  unnest(scores) %>%
  separate(fraction, into = c("control", "fraction"), sep = "%%") %>%
  mutate(fraction = as.numeric(fraction))

# Filter signatures
sigs_filt3 <- scores_all_sims_tidy %>%
  group_by(celltype, signature, sim_id) %>%
  summarise(cor = cor(fraction, score, method = "pearson")) %>%
  group_by(celltype, signature) %>%
  summarise(mean_cor = mean(cor)) %>%
  top_frac(n=0.1, wt=mean_cor) %>%
  pull(signature)

