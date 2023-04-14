###############################################################################
# This script run xCell and xCell2 and load Kassandra's benchmarking results
###############################################################################

# Input - xCell2.0 signatures RDS
# xcell2_sigs <- xCell2Ref.S4
# saveRDS(xcell2_sigs, "/home/almogangel/xCell2_git/Data/benchmarking_data/xcell2ref_sref_main.rds")
xcell2_sigs <- readRDS("/home/almogangel/xCell2_git/Data/benchmarking_data/xcell2ref_ts_main.rds")
outfile <- "/home/almogangel/xCell2_git/Data/benchmarking_data/kass_benchmarking_cors_ts_main.rds"

library(tidyverse)
source("/home/almogangel/xCell2_git/R/xCell2.R")
xCell.data <- xCell::xCell.data


# Load celltype conversion
celltype_conversion_long <- read_tsv("/home/almogangel/xCell2_git/Data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))


# Load benchmarking truths and mixtures
truths_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/"
mix_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/"
ds <- gsub(".tsv", "", list.files(truths_dir))


# This function returns the correlation between each method results to the ground truth
getCorrelations <- function(dataset){

  # Load mixture
  mix <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/", dataset, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE)

  # Load truth
  truth <- read_tsv(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/", dataset, ".tsv")) %>%
    dplyr::rename(celltype = 1) %>% # Fix for some datasets
    mutate(celltype = plyr::mapvalues(celltype, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE))

  samples2use <- intersect(colnames(truth)[-1], colnames(mix))
  mix <- mix[,samples2use]
  truth <- truth[,c(colnames(truth)[1], samples2use)]

  if(!all(colnames(mix) == colnames(truth)[-1])){
    print(paste0("Problem with ds:" , file))
    break
  }

  truth <- pivot_longer(truth, -celltype, names_to = "sample", values_to = "true_fracs")

  # Run xCell2
  xcell2.out <- xCell2Analysis(mix, xcell2ref = xcell2_sigs)

  # Run xCell
  xcell.out <- as.data.frame(xCell::xCellAnalysis(mix,rnaseq = FALSE))
  rownames(xcell.out) <- plyr::mapvalues(rownames(xcell.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

  # Load results from other methods
  kass.out <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/Kassandra/", dataset ,"_predicted_by_Kassandra.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(kass.out) <- plyr::mapvalues(rownames(kass.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  abis.out <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/ABIS/", dataset, "_predicted_by_ABIS.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(abis.out) <- plyr::mapvalues(rownames(abis.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  epic_bref.out <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/EPIC_BRef/", dataset, "_predicted_by_EPIC_BRef.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(epic_bref.out) <- plyr::mapvalues(rownames(epic_bref.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  epic_tref.out <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/EPIC_TRef/", dataset, "_predicted_by_EPIC_TRef.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(epic_tref.out) <- plyr::mapvalues(rownames(epic_tref.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  fardeep_abs.in <- paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/fardeep_absolute/", dataset, "_predicted_by_fardeep_absolute.tsv")
  if(file.exists(fardeep_abs.in)){
    fardeep_abs.out <- read.table(fardeep_abs.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(fardeep_abs.out) <- plyr::mapvalues(rownames(fardeep_abs.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    fardeep_abs.out <- NULL
  }
  fardeep_rel.in <- paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/fardeep_relative/", dataset, "_predicted_by_fardeep_relative.tsv")
  if(file.exists(fardeep_rel.in)){
    fardeep_rel.out <- read.table(fardeep_rel.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(fardeep_rel.out) <- plyr::mapvalues(rownames(fardeep_rel.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    fardeep_rel.out <- NULL
  }
  cbrx_hnsc.out <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/HNSC_matrix_with_BC_CIBERSORTX/", dataset, "_predicted_by_HNSC_matrix_with_BC_CIBERSORTX.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(cbrx_hnsc.out) <- plyr::mapvalues(rownames(cbrx_hnsc.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  cbrx_lm22.out <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/LM22_matrix_with_BC_CIBERSORTX/", dataset, "_predicted_by_LM22_matrix_with_BC_CIBERSORTX.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  cbrx_lm22.out[rownames(cbrx_lm22.out) == "Dendritic cells resting",] <- colSums(cbrx_lm22.out[rownames(cbrx_lm22.out) %in% c("Dendritic cells resting", "Dendritic cells activated"),])
  rownames(cbrx_lm22.out[rownames(cbrx_lm22.out) == "Dendritic cells resting",]) <- "DC"
  cbrx_lm22.out <- cbrx_lm22.out[!rownames(cbrx_lm22.out) %in% c("Mast cells resting", "Mast cells activated", "NK cells resting", "NK cells activated",
                                                                 "Dendritic cells resting", "Dendritic cells activated"),] # Fix for some datasets
  rownames(cbrx_lm22.out) <- plyr::mapvalues(rownames(cbrx_lm22.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  quanti.out <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/quantiseq/", dataset, "_predicted_by_quantiseq.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(quanti.out) <- plyr::mapvalues(rownames(quanti.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  quanti_tumor.out <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/quantiseq_tumor/", dataset, "_predicted_by_quantiseq_tumor.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(quanti_tumor.out) <- plyr::mapvalues(rownames(quanti_tumor.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

  scaden.in <- paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/Scaden_PBMC/", dataset, "_predicted_by_Scaden_PBMC.tsv")
  if(file.exists(scaden.in)){
    scaden.out <- read.table(scaden.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(scaden.out) <- plyr::mapvalues(rownames(scaden.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    scaden.out <- NULL
  }


  all_methods <- c("xCell2.0", "xCell", "Kassandara", "ABIS", "EPIC-BRef", "EPIC-TRef", "FARDEEP-Abs", "FARDEEP-Rel",
                   "CIBERSORTx-HNSC", "CIBERSORTx-LM22", "quanTIseq", "quanTIseq-T", "Scaden")

  all.out <- list(xcell2.out, xcell.out, kass.out, abis.out, epic_bref.out, epic_tref.out, fardeep_abs.out,
                  fardeep_rel.out, cbrx_hnsc.out, cbrx_lm22.out, quanti.out, quanti_tumor.out, scaden.out)
  names(all.out) <- all_methods

  all.out <- compact(all.out)
  shared_samples <- Reduce(intersect, lapply(all.out, names))
  # shared_celltypes <- Reduce(intersect, lapply(all.out, rownames))
  # sort(table(unlist(lapply(all.out, rownames))), decreasing = T)
  # shared_celltypes <- unique(c(shared_celltypes, celltypes2add))
  #shared_celltypes <- unique(unlist(lapply(all.out, rownames)))


  all.out <- lapply(all.out, function(x){x[, shared_samples]})
  all.out <- lapply(all.out, function(x){cbind(celltype = rownames(x), x)})
  all.out <- lapply(all.out, function(x){drop_na(x)})
  all.out <- lapply(all.out, function(x){x[rowSums(x != 0) > 3,]}) # Minimum 3 values for cell type

  cors.out <- enframe(all.out, name = "method", value = "predictions") %>%
    unnest(cols = c(predictions)) %>%
    pivot_longer(-c(method, celltype), names_to = "sample", values_to = "prediction") %>%
    left_join(truth, by = c("celltype", "sample")) %>%
    mutate(true_fracs = as.numeric(true_fracs)) %>%
    drop_na() %>%
    group_by(method, celltype) %>%
    summarise(Pearson = cor(prediction, true_fracs),
              Spearman = cor(prediction, true_fracs, method = "spearman"),
              MSE = mean((prediction - true_fracs)^2)) %>%
    pivot_longer(cols = c("Pearson", "Spearman"), names_to = "cor_method", values_to = "cor_score") %>%
    mutate(dataset = dataset)

  return(cors.out)
}


# Use benchmarking for those datasets
cors.out.list <- vector("list", length(ds))
names(cors.out.list) <- ds


for (d in ds) {

  # print(paste0(d, " --------------------------------------------------------- "))
  # cors.out <- getCorrelations(dataset = d) %>%
  #   mutate(dataset = d)
  # cors.out.list[[d]] <- cors.out

  print(paste0(d, " --------------------------------------------------------- "))
  cors.out.list[[d]] <- getCorrelations(dataset = d)


}


saveRDS(cors.out.list, outfile)
