############################################################
# Using only protein coding genes helps?
############################################################


library(tidyverse)
library(xCell2)
library(parallel)

setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")


# "/bigdata/almogangel/xCell2/dev_scripts/prep_ref_val_pairs.R"
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
sc.refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val.rds")
# Load validation data
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
sc.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc.vals.rds")

refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))


# Load references matrices
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


vals.refs.res <- refval.tbl %>%
  # Get number of samples in the validation dataset
  mutate(n_val_samples = ncol(cyto.vals$truth[[val_type]][[val_dataset]])) %>%
  filter(n_shared_celltypes > 2) %>%
  mutate(method = "xCell2", .before = everything())


vals.refs.res.sc <- sc.refval.tbl %>%
  # Get number of samples in the validation dataset
  mutate(n_val_samples = ncol(sc.vals$truth[[val_type]][[val_dataset]])) %>%
  filter(n_shared_celltypes > 2) %>%
  mutate(method = "xCell2", .before = everything())



# Run settings
set.seed(123)
thisseed <- 123
cores2use <- 20

# gene settings
useTopVar <- TRUE
nTopVar <- 5000
#ProteinCoding <- FALSE
#ProteinCodingSC <- TRUE
genesGroups <- c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY")
genesGroupsSC <- c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY")


# signature settings
load_sigs <- FALSE
sigs_suffix <- "x"
scores_results <- TRUE
minpbcells <- 30
minpbgroups <- 10
weight_genes <- TRUE

# simulations settings
simNoise <- 0
gamma2use <- 0.5

# xCell2Analysis
tranform <- TRUE
spillover <- FALSE


# Cytometry validations ----------------------
ProteinCoding <- TRUE
ProteinCodingSC <- TRUE
bulkTscT_res <- mclapply(1:nrow(vals.refs.res), function(i) {

  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }


  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }



  if (scores_results) {
    # Return scores results
    if (is.null(sigsFile)) {
      sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                          sigsFile = NULL, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    }else{
      sigs <- readRDS(sigsFile)
    }

    mix_ranked <- singscore::rankGenes(mix)
    sigs.scores <- sapply(sigs, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    sigs.scores <- t(sigs.scores)
    colnames(sigs.scores) <- colnames(mix_ranked)
    sigs.scores

    celltypes <- vals.refs.res[i,]$shared_celltypes[[1]]

    t(sapply(celltypes, function(ct){
      ct_index <- gsub("#.*", "", rownames(sigs.scores)) == ct
      colMeans(sigs.scores[ct_index,])
    }))

  }else{
    # Return model results

    # xCell2Train
    sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                                 sim_noise = simNoise, regGamma = gamma2use, nCores = cores2use)
    print("sigs object ready.")

    # xCell2Analysis
    res <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
    res <- res[vals.refs.res[i,]$shared_celltypes[[1]], ]
    res

  }


}, mc.cores = 5, mc.set.seed = FALSE)

ProteinCoding <- FALSE
ProteinCodingSC <- TRUE
bulkFscT_res <- mclapply(1:nrow(vals.refs.res), function(i) {

  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }


  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }



  if (scores_results) {
    # Return scores results
    if (is.null(sigsFile)) {
      sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                          sigsFile = NULL, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    }else{
      sigs <- readRDS(sigsFile)
    }

    mix_ranked <- singscore::rankGenes(mix)
    sigs.scores <- sapply(sigs, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    sigs.scores <- t(sigs.scores)
    colnames(sigs.scores) <- colnames(mix_ranked)
    sigs.scores

    celltypes <- vals.refs.res[i,]$shared_celltypes[[1]]

    t(sapply(celltypes, function(ct){
      ct_index <- gsub("#.*", "", rownames(sigs.scores)) == ct
      colMeans(sigs.scores[ct_index,])
    }))

  }else{
    # Return model results

    # xCell2Train
    sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                                 sim_noise = simNoise, regGamma = gamma2use, nCores = cores2use)
    print("sigs object ready.")

    # xCell2Analysis
    res <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
    res <- res[vals.refs.res[i,]$shared_celltypes[[1]], ]
    res

  }


}, mc.cores = 5, mc.set.seed = FALSE)

ProteinCoding <- FALSE
ProteinCodingSC <- FALSE
bulkFscF_res <- mclapply(1:nrow(vals.refs.res), function(i) {

  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }


  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }



  if (scores_results) {
    # Return scores results
    if (is.null(sigsFile)) {
      sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                          sigsFile = NULL, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    }else{
      sigs <- readRDS(sigsFile)
    }

    mix_ranked <- singscore::rankGenes(mix)
    sigs.scores <- sapply(sigs, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    sigs.scores <- t(sigs.scores)
    colnames(sigs.scores) <- colnames(mix_ranked)
    sigs.scores

    celltypes <- vals.refs.res[i,]$shared_celltypes[[1]]

    t(sapply(celltypes, function(ct){
      ct_index <- gsub("#.*", "", rownames(sigs.scores)) == ct
      colMeans(sigs.scores[ct_index,])
    }))

  }else{
    # Return model results

    # xCell2Train
    sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                                 sim_noise = simNoise, regGamma = gamma2use, nCores = cores2use)
    print("sigs object ready.")

    # xCell2Analysis
    res <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
    res <- res[vals.refs.res[i,]$shared_celltypes[[1]], ]
    res

  }


}, mc.cores = 5, mc.set.seed = FALSE)

ProteinCoding <- TRUE
ProteinCodingSC <- FALSE
bulkTscF_res <- mclapply(1:nrow(vals.refs.res), function(i) {

  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }


  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }



  if (scores_results) {
    # Return scores results
    if (is.null(sigsFile)) {
      sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                          sigsFile = NULL, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    }else{
      sigs <- readRDS(sigsFile)
    }

    mix_ranked <- singscore::rankGenes(mix)
    sigs.scores <- sapply(sigs, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    sigs.scores <- t(sigs.scores)
    colnames(sigs.scores) <- colnames(mix_ranked)
    sigs.scores

    celltypes <- vals.refs.res[i,]$shared_celltypes[[1]]

    t(sapply(celltypes, function(ct){
      ct_index <- gsub("#.*", "", rownames(sigs.scores)) == ct
      colMeans(sigs.scores[ct_index,])
    }))

  }else{
    # Return model results

    # xCell2Train
    sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                                 sim_noise = simNoise, regGamma = gamma2use, nCores = cores2use)
    print("sigs object ready.")

    # xCell2Analysis
    res <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
    res <- res[vals.refs.res[i,]$shared_celltypes[[1]], ]
    res

  }


}, mc.cores = 5, mc.set.seed = FALSE)


getCors <- function(res_mat, valList, valType, valDataset, corMethod = "spearman"){

  truth_mat <- valList$truth[[valType]][[valDataset]]


  celltypes <- intersect(rownames(res_mat), rownames(truth_mat))

  cor.tests <- lapply(celltypes, function(ct){

    truth <- truth_mat[ct,]

    if(any(is.na(truth))){
      stop("NAs in truth!!!!")
    }

    res <- res_mat[ct,]

    samples <- intersect(names(res), names(truth))

    cor.res <- suppressWarnings(cor(res[samples], truth[samples], method = "spearman"))
    cor.res <- ifelse(is.na(cor.res), 0, cor.res)
    cor.res

  })
  names(cor.tests) <- celltypes

  return(cor.tests)
}


# <param> <reference> <validation> <n_samples> <celltype> <cor>
bulkTscT_res.cors <- vals.refs.res %>%
  ungroup() %>%
  select(method, val_type, val_dataset, n_val_samples) %>%
  mutate(method = "only protein coding",
         res = bulkTscT_res) %>%
  rowwise() %>%
  mutate(cors = list(getCors(res, valList = cyto.vals, valDataset = val_dataset, valType = val_type))) %>%
  select(-res) %>%
  unnest_longer(cors, values_to = "cor", indices_to = "celltype")

bulkFscT_res.cors <- vals.refs.res %>%
  ungroup() %>%
  select(method, val_type, val_dataset, n_val_samples) %>%
  mutate(method = "protein coding for SC",
         res = bulkFscT_res) %>%
  rowwise() %>%
  mutate(cors = list(getCors(res, valList = cyto.vals, valDataset = val_dataset, valType = val_type))) %>%
  select(-res) %>%
  unnest_longer(cors, values_to = "cor", indices_to = "celltype")

bulkFscF_res.cors <- vals.refs.res %>%
  ungroup() %>%
  select(method, val_type, val_dataset, n_val_samples) %>%
  mutate(method = "no filtering for protein coding",
         res = bulkFscF_res) %>%
  rowwise() %>%
  mutate(cors = list(getCors(res, valList = cyto.vals, valDataset = val_dataset, valType = val_type))) %>%
  select(-res) %>%
  unnest_longer(cors, values_to = "cor", indices_to = "celltype")

bulkTscF_res.cors <- vals.refs.res %>%
  ungroup() %>%
  select(method, val_type, val_dataset, n_val_samples) %>%
  mutate(method = "protein coding for bulk",
         res = bulkTscF_res) %>%
  rowwise() %>%
  mutate(cors = list(getCors(res, valList = cyto.vals, valDataset = val_dataset, valType = val_type))) %>%
  select(-res) %>%
  unnest_longer(cors, values_to = "cor", indices_to = "celltype")


allcors <- rbind(bulkTscT_res.cors, bulkFscT_res.cors, bulkFscF_res.cors, bulkTscF_res.cors)

saveRDS(allcors, "/bigdata/almogangel/xCell2/dev_scripts/testing_params/protein_coding_genes_allcors_cyto.rds")
allcors <- readRDS("/bigdata/almogangel/xCell2/dev_scripts/testing_params/protein_coding_genes_allcors_cyto.rds")

allcors <- allcors %>%
  group_by(method) %>%
  mutate(median_cor = median(cor, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(median_cor))

ggplot(allcors, aes(x = reorder(method, -median_cor), y = cor, fill=method)) +
  geom_boxplot() +
  geom_jitter(width = 0.3, alpha = 0.6) + # adding jitter
  labs(x = "", y = "Spearman Rho", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # rotating x-axis labels for better readability



# Single cell validations ----------------------
ProteinCoding <- TRUE
ProteinCodingSC <- TRUE
bulkTscT_res.sc <- mclapply(1:nrow(vals.refs.res.sc), function(i) {

  # Load data
  val_ref <- paste0(vals.refs.res.sc[i,]$val_dataset, "_", vals.refs.res.sc[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- sc.vals$mixtures[[vals.refs.res.sc[i,]$val_type]][[vals.refs.res.sc[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res.sc[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res.sc[i,]$ref_type)


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }


  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }



  if (scores_results) {
    # Return scores results
    if (is.null(sigsFile)) {
      sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                          sigsFile = NULL, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    }else{
      sigs <- readRDS(sigsFile)
    }

    mix_ranked <- singscore::rankGenes(mix)
    sigs.scores <- sapply(sigs, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    sigs.scores <- t(sigs.scores)
    colnames(sigs.scores) <- colnames(mix_ranked)
    sigs.scores

    celltypes <- vals.refs.res.sc[i,]$shared_celltypes[[1]]

    t(sapply(celltypes, function(ct){
      ct_index <- gsub("#.*", "", rownames(sigs.scores)) == ct
      colMeans(sigs.scores[ct_index,])
    }))

  }else{
    # Return model results

    # xCell2Train
    sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                                 sim_noise = simNoise, regGamma = gamma2use, nCores = cores2use)
    print("sigs object ready.")

    # xCell2Analysis
    res <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
    res <- res[vals.refs.res.sc[i,]$shared_celltypes[[1]], ]
    res

  }


}, mc.cores = 5, mc.set.seed = FALSE)

ProteinCoding <- FALSE
ProteinCodingSC <- TRUE
bulkFscT_res.sc <- mclapply(1:nrow(vals.refs.res.sc), function(i) {

  # Load data
  val_ref <- paste0(vals.refs.res.sc[i,]$val_dataset, "_", vals.refs.res.sc[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- sc.vals$mixtures[[vals.refs.res.sc[i,]$val_type]][[vals.refs.res.sc[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res.sc[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res.sc[i,]$ref_type)


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }


  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }



  if (scores_results) {
    # Return scores results
    if (is.null(sigsFile)) {
      sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                          sigsFile = NULL, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    }else{
      sigs <- readRDS(sigsFile)
    }

    mix_ranked <- singscore::rankGenes(mix)
    sigs.scores <- sapply(sigs, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    sigs.scores <- t(sigs.scores)
    colnames(sigs.scores) <- colnames(mix_ranked)
    sigs.scores

    celltypes <- vals.refs.res.sc[i,]$shared_celltypes[[1]]

    t(sapply(celltypes, function(ct){
      ct_index <- gsub("#.*", "", rownames(sigs.scores)) == ct
      colMeans(sigs.scores[ct_index,])
    }))

  }else{
    # Return model results

    # xCell2Train
    sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                                 sim_noise = simNoise, regGamma = gamma2use, nCores = cores2use)
    print("sigs object ready.")

    # xCell2Analysis
    res <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
    res <- res[vals.refs.res.sc[i,]$shared_celltypes[[1]], ]
    res

  }


}, mc.cores = 5, mc.set.seed = FALSE)

ProteinCoding <- FALSE
ProteinCodingSC <- FALSE
bulkFscF_res.sc <- mclapply(1:nrow(vals.refs.res.sc), function(i) {

  # Load data
  val_ref <- paste0(vals.refs.res.sc[i,]$val_dataset, "_", vals.refs.res.sc[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- sc.vals$mixtures[[vals.refs.res.sc[i,]$val_type]][[vals.refs.res.sc[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res.sc[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res.sc[i,]$ref_type)


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }


  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }



  if (scores_results) {
    # Return scores results
    if (is.null(sigsFile)) {
      sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                          sigsFile = NULL, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    }else{
      sigs <- readRDS(sigsFile)
    }

    mix_ranked <- singscore::rankGenes(mix)
    sigs.scores <- sapply(sigs, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    sigs.scores <- t(sigs.scores)
    colnames(sigs.scores) <- colnames(mix_ranked)
    sigs.scores

    celltypes <- vals.refs.res.sc[i,]$shared_celltypes[[1]]

    t(sapply(celltypes, function(ct){
      ct_index <- gsub("#.*", "", rownames(sigs.scores)) == ct
      colMeans(sigs.scores[ct_index,])
    }))

  }else{
    # Return model results

    # xCell2Train
    sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                                 sim_noise = simNoise, regGamma = gamma2use, nCores = cores2use)
    print("sigs object ready.")

    # xCell2Analysis
    res <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
    res <- res[vals.refs.res.sc[i,]$shared_celltypes[[1]], ]
    res

  }


}, mc.cores = 5, mc.set.seed = FALSE)

ProteinCoding <- TRUE
ProteinCodingSC <- FALSE
bulkTscF_res.sc <- mclapply(1:nrow(vals.refs.res.sc), function(i) {

  # Load data
  val_ref <- paste0(vals.refs.res.sc[i,]$val_dataset, "_", vals.refs.res.sc[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- sc.vals$mixtures[[vals.refs.res.sc[i,]$val_type]][[vals.refs.res.sc[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res.sc[i,]$ref_type]][[vals.refs.res.sc[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res.sc[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res.sc[i,]$ref_type)


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }


  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }



  if (scores_results) {
    # Return scores results
    if (is.null(sigsFile)) {
      sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                          sigsFile = NULL, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    }else{
      sigs <- readRDS(sigsFile)
    }

    mix_ranked <- singscore::rankGenes(mix)
    sigs.scores <- sapply(sigs, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    sigs.scores <- t(sigs.scores)
    colnames(sigs.scores) <- colnames(mix_ranked)
    sigs.scores

    celltypes <- vals.refs.res.sc[i,]$shared_celltypes[[1]]

    t(sapply(celltypes, function(ct){
      ct_index <- gsub("#.*", "", rownames(sigs.scores)) == ct
      colMeans(sigs.scores[ct_index,])
    }))

  }else{
    # Return model results

    # xCell2Train
    sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                                 sim_noise = simNoise, regGamma = gamma2use, nCores = cores2use)
    print("sigs object ready.")

    # xCell2Analysis
    res <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
    res <- res[vals.refs.res.sc[i,]$shared_celltypes[[1]], ]
    res

  }


}, mc.cores = 5, mc.set.seed = FALSE)


getCors <- function(res_mat, valList, valType, valDataset, corMethod = "spearman"){

  truth_mat <- valList$truth[[valType]][[valDataset]]


  celltypes <- intersect(rownames(res_mat), rownames(truth_mat))

  cor.tests <- lapply(celltypes, function(ct){

    truth <- truth_mat[ct,]

    if(any(is.na(truth))){
      stop("NAs in truth!!!!")
    }

    res <- res_mat[ct,]

    samples <- intersect(names(res), names(truth))

    cor.res <- suppressWarnings(cor(res[samples], truth[samples], method = "spearman"))
    cor.res <- ifelse(is.na(cor.res), 0, cor.res)
    cor.res

  })
  names(cor.tests) <- celltypes

  return(cor.tests)
}


# <param> <reference> <validation> <n_samples> <celltype> <cor>
bulkTscT_res.cors.sc <- vals.refs.res.sc %>%
  ungroup() %>%
  select(method, val_type, val_dataset, n_val_samples) %>%
  mutate(method = "only protein coding",
         res = bulkTscT_res.sc) %>%
  rowwise() %>%
  mutate(cors = list(getCors(res, valList = sc.vals, valDataset = val_dataset, valType = val_type))) %>%
  select(-res) %>%
  unnest_longer(cors, values_to = "cor", indices_to = "celltype")

bulkFscT_res.cors.sc <- vals.refs.res.sc %>%
  ungroup() %>%
  select(method, val_type, val_dataset, n_val_samples) %>%
  mutate(method = "protein coding for SC",
         res = bulkFscT_res.sc) %>%
  rowwise() %>%
  mutate(cors = list(getCors(res, valList = sc.vals, valDataset = val_dataset, valType = val_type))) %>%
  select(-res) %>%
  unnest_longer(cors, values_to = "cor", indices_to = "celltype")

bulkFscF_res.cors.sc <- vals.refs.res.sc %>%
  ungroup() %>%
  select(method, val_type, val_dataset, n_val_samples) %>%
  mutate(method = "no filtering for protein coding",
         res = bulkFscF_res.sc) %>%
  rowwise() %>%
  mutate(cors = list(getCors(res, valList = sc.vals, valDataset = val_dataset, valType = val_type))) %>%
  select(-res) %>%
  unnest_longer(cors, values_to = "cor", indices_to = "celltype")

bulkTscF_res.cors.sc <- vals.refs.res.sc %>%
  ungroup() %>%
  select(method, val_type, val_dataset, n_val_samples) %>%
  mutate(method = "protein coding for bulk",
         res = bulkTscF_res.sc) %>%
  rowwise() %>%
  mutate(cors = list(getCors(res, valList = sc.vals, valDataset = val_dataset, valType = val_type))) %>%
  select(-res) %>%
  unnest_longer(cors, values_to = "cor", indices_to = "celltype")


allcors <- rbind(bulkTscT_res.cors.sc, bulkFscT_res.cors.sc, bulkFscF_res.cors.sc, bulkTscF_res.cors.sc)


saveRDS(allcors, "/bigdata/almogangel/xCell2/dev_scripts/testing_params/protein_coding_genes_allcors_sc.rds")


# allcors <- allcors %>%
#   group_by(method) %>%
#   mutate(median_cor = median(cor, na.rm = TRUE)) %>%
#   ungroup() %>%
#   arrange(desc(median_cor))
#
# ggplot(allcors, aes(x = reorder(method, -median_cor), y = cor, fill=method)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.3, alpha = 0.6) + # adding jitter
#   labs(x = "", y = "Spearman Rho", fill = "") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) # rotating x-axis labels for better readability
#
