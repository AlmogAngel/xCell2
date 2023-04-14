###############################################################################
# This script run xCell, xCell2 and CIBERSORTx with the same reference
###############################################################################

library(tidyverse)

source("R/xCell2.R")
source("dev_scripts/run_cibersortx.R")
xCell.data <- xCell::xCell.data
celltype_conversion <- read_tsv("Data/celltype_conversion_with_ontology.txt")
celltype_conversion_long <- celltype_conversion %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

getCorrelation <- function(results_mat, truth, cor_method = "spearman"){

  samples2use <- intersect(colnames(truth), colnames(results_mat))
  results_mat <- results_mat[,samples2use]
  truth <- truth[,samples2use]

  if (!all(colnames(truth) == colnames(results_mat))) {
    errorCondition(paste0("Error with dataset: ", file))
  }

  celltypes2use <- intersect(rownames(truth), rownames(results_mat))

  all_celltypes_cor <- sapply(celltypes2use, function(ctoi){

    scores_ctoi <- results_mat[ctoi,]
    truth_ctoi <- truth[ctoi, names(scores_ctoi)]
    truth_ctoi <- truth_ctoi[which(!is.na(truth_ctoi))]
    truth_ctoi <- truth_ctoi[which(truth_ctoi != "")]
    scores_ctoi <- scores_ctoi[names(truth_ctoi)]

    if (!all(names(truth_ctoi) == colnames(scores_ctoi))) {
      errorCondition(paste0("Error with dataset: ", file))
    }

    if (all(as.numeric(truth_ctoi) == 0) | all(as.numeric(scores_ctoi) == 0)) {
      NULL
    }else{
      cor(as.numeric(scores_ctoi), as.numeric(truth_ctoi), method = cor_method)
    }
  })

  if (is.list(all_celltypes_cor)) {
    all_celltypes_cor <- unlist(all_celltypes_cor)
  }

  out <- enframe(compact(all_celltypes_cor)) %>%
    arrange(-value) %>%
    rename(celltype = name , cor = value)

  return(out)
}


# xcell2_blood_sigs <- xCell2Train(ref = ref, labels = labels, data_type = "rnaseq", ontology_file_checked = ontology_file_checked)
# saveRDS(xcell2_blood_sigs, "/Users/almogang/Documents/xCell2.0/reference_data/xcell2_blood_sigs_simulations.rds")

# Blood validation ----

# Single cell reference (Tabula Sapiens) ----
xcell2_blood_sigs <- readRDS("/Users/almogang/Documents/xCell2.0/reference_data/xcell2ref_ts_blood_main.rds")
lm22 <- "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/LM22.txt"

# Load blood validation datasets
blood_ds <- c("BG_blood", "GSE107011", "GSE107572", "GSE115823", "GSE127813", "GSE53655", "GSE60424")


# Run benchmarking
blood_validation.list <- list()
for (ds in blood_ds) {

  print(paste0("Working on: ", ds))

  # Load mixture and true fractions
  mix_path <- paste0("/Users/almogang/Documents/xCell2.0/validation_data/expressions/", ds, "_expressions.tsv")
  mix <- read.table(mix_path, check.names = FALSE, row.names = 1, header = TRUE)
  truth <- read.table(paste0("/Users/almogang/Documents/xCell2.0/validation_data/cell_values/", ds, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

  # Run xCell2.0
  xcell2.out <- xCell2Analysis(mix = mix, xcell2ref = xcell2_blood_sigs)
  print("xCell2.0 - Done")

  # Run xCell
  xcell.out <- as.data.frame(xCell::xCellAnalysis(mix, rnaseq = TRUE))
  rownames(xcell.out) <- plyr::mapvalues(rownames(xcell.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  print("xCell - Done")

  # Run CIBERSORTx - custom reference
  cbrx_custom.out <- runCIBERSORTxSC(mix = mix_path, ref = cbrx_blood_sig_mat, absolute = FALSE, rmbatchBmode = TRUE)
  print("CIBERSORTx (custom reference) - Done")

  # Run CIBERSORTx - LM22
  cbrx_lm22.out <- runCIBERSORTx(mix = mix_path, sig_mat = lm22, absolute = FALSE, rmbatchBmode = TRUE)
  # Add broad cell types to LM22
  rownames(cbrx_lm22.out) <- plyr::mapvalues(rownames(cbrx_lm22.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  tmp <- t(data.frame("B-cells" = colSums(cbrx_lm22.out[c("naive B-cells", "Memory B-cells", "Plasma cells"),]),
                      "CD4+ T-cells" = colSums(cbrx_lm22.out[c("CD4+ naive T-cells", "T cells CD4 memory resting", "CD4+ memory T-cells"),]),
                      "Macrophages" = colSums(cbrx_lm22.out[c("Macrophages M0", "Macrophages M1", "Macrophages M2"),]),
                      "T-cells" = colSums(cbrx_lm22.out[c("CD8+ T-cells", "CD4+ naive T-cells", "T cells CD4 memory resting", "CD4+ memory T-cells",
                                                          "Follicular T-helper", "Tregs", "Tgd cells"),]), check.names = FALSE))
  cbrx_lm22.out <- rbind(cbrx_lm22.out, tmp)
  print("CIBERSORTx (LM22) - Done")


  # Calculate Correlations
  xcell2.cor <- getCorrelation(results_mat = xcell2.out, truth = truth, cor_method = "spearman") %>%
    mutate(method = "xCell2.0 - Kassandra's Blood")
  xcell.cor <- getCorrelation(results_mat = xcell.out, truth = truth, cor_method = "spearman") %>%
    mutate(method = "xCell")
  cbrx_custom.cor <- getCorrelation(results_mat = cbrx_custom.out, truth = truth, cor_method = "spearman") %>%
    mutate(method = "CIBERSORTx - Kassandra's Blood")
  cbrx_lm22.cor <- getCorrelation(results_mat = cbrx_lm22.out, truth = truth, cor_method = "spearman") %>%
    mutate(method = "CIBERSORTx - LM22")

  blood_validation.list[[ds]] <- rbind(xcell2.cor, xcell.cor, cbrx_custom.cor, cbrx_lm22.cor)

}


# Bulk human reference  ----
# Load signatures for xCell2.0 and signatures matrix for CIBERSORTx
xcell2_blood_sigs <- readRDS("/Users/almogang/Documents/xCell2.0/reference_data/xcell2_blood_sigs_simulations.rds")
cbrx_blood_sig_mat <- "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/cibersortx_blood_ref.txt-sigmat/CIBERSORTx_cibersortx_blood_pheno_df.CIBERSORTx_cibersortx_blood_ref.bm.K999.txt"
lm22 <- "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/LM22.txt"

# Load blood validation datasets
blood_ds <- c("BG_blood", "GSE107011", "GSE107572", "GSE115823", "GSE127813", "GSE53655", "GSE60424")
# "GSE64655", "sc_pbmc", "SDY67"

# Run benchmarking
blood_validation.list <- list()
for (ds in blood_ds) {

  print(paste0("Working on: ", ds))

  # Load mixture and true fractions
  mix_path <- paste0("/Users/almogang/Documents/xCell2.0/validation_data/expressions/", ds, "_expressions.tsv")
  mix <- read.table(mix_path, check.names = FALSE, row.names = 1, header = TRUE)
  truth <- read.table(paste0("/Users/almogang/Documents/xCell2.0/validation_data/cell_values/", ds, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

  # Run xCell2.0
  xcell2.out <- xCell2Analysis(mix = mix, xcell2ref = xcell2_blood_sigs)
  print("xCell2.0 - Done")

  # Run xCell
  xcell.out <- as.data.frame(xCell::xCellAnalysis(mix, rnaseq = TRUE))
  rownames(xcell.out) <- plyr::mapvalues(rownames(xcell.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  print("xCell - Done")

  # Run CIBERSORTx - custom reference
  cbrx_custom.out <- runCIBERSORTx(mix = mix_path, sig_mat = cbrx_blood_sig_mat, absolute = FALSE, rmbatchBmode = TRUE)
  print("CIBERSORTx (custom reference) - Done")

  # Run CIBERSORTx - LM22
  cbrx_lm22.out <- runCIBERSORTx(mix = mix_path, sig_mat = lm22, absolute = FALSE, rmbatchBmode = TRUE)
  # Add broad cell types to LM22
  rownames(cbrx_lm22.out) <- plyr::mapvalues(rownames(cbrx_lm22.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  tmp <- t(data.frame("B-cells" = colSums(cbrx_lm22.out[c("naive B-cells", "Memory B-cells", "Plasma cells"),]),
            "CD4+ T-cells" = colSums(cbrx_lm22.out[c("CD4+ naive T-cells", "T cells CD4 memory resting", "CD4+ memory T-cells"),]),
            "Macrophages" = colSums(cbrx_lm22.out[c("Macrophages M0", "Macrophages M1", "Macrophages M2"),]),
            "T-cells" = colSums(cbrx_lm22.out[c("CD8+ T-cells", "CD4+ naive T-cells", "T cells CD4 memory resting", "CD4+ memory T-cells",
                                                        "Follicular T-helper", "Tregs", "Tgd cells"),]), check.names = FALSE))
  cbrx_lm22.out <- rbind(cbrx_lm22.out, tmp)
  print("CIBERSORTx (LM22) - Done")


  # Calculate Correlations
  xcell2.cor <- getCorrelation(results_mat = xcell2.out, truth = truth, cor_method = "spearman") %>%
    mutate(method = "xCell2.0 - Kassandra's Blood")
  xcell.cor <- getCorrelation(results_mat = xcell.out, truth = truth, cor_method = "spearman") %>%
    mutate(method = "xCell")
  cbrx_custom.cor <- getCorrelation(results_mat = cbrx_custom.out, truth = truth, cor_method = "spearman") %>%
    mutate(method = "CIBERSORTx - Kassandra's Blood")
  cbrx_lm22.cor <- getCorrelation(results_mat = cbrx_lm22.out, truth = truth, cor_method = "spearman") %>%
    mutate(method = "CIBERSORTx - LM22")

  blood_validation.list[[ds]] <- rbind(xcell2.cor, xcell.cor, cbrx_custom.cor, cbrx_lm22.cor)

}

# Plot ----

lapply(blood_validation.list, function(x) unique(pull(x, celltype)))
blood_validation.list.sub <- blood_validation.list[!names(blood_validation.list) %in% c("GSE53655", "GSE60424", "GSE115823")]


blood_validation.tbl <- enframe(blood_validation.list.sub) %>%
  unnest(cols = c(value)) %>%
  rename(dataset = name)


ct2plot <- c("Monocytes", "Neutrophils", "NK cells", "cDC", "T-cells", "Tregs",
             "Memory T-helpers", "CD4+ T-cells", "CD8+ T-cells", "B-cells",
             "Fibroblasts", "Lymphocytes", "Cancer cells", "CD8+ T-cells PD1 low", "CD8+ T-cells PD1 high",
             "Macrophages", "Macrophages M1", "Macrophages M2", "Endothelial cells", "Non plasma B-cells")

ct2plot <- Reduce(intersect, lapply(blood_validation.list.sub, function(x) pull(x, celltype)))


blood_validation.tbl %>%
  filter(celltype %in% ct2plot) %>%
  mutate(is_xcell2 = ifelse(startsWith(method, "xCell2.0"), "yes", "no")) %>%
  ggplot(., aes(x=method, y= cor)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=celltype), size=4, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10,12)],
                              "#424242", "#8B1C62", "#00F5FF", "#FF3E96")) +
                                scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 10, angle = 30, hjust=1),
        axis.text.y = element_text(size = 12, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r (median)", title = "BG_blood Validation Dataset", x = NULL, colour = NULL, fill = NULL) +
  facet_wrap(~dataset, scales = "free")





# old -------------------------------------------

blood_ds <- c("BG_blood", "GSE107011", "GSE107572", "GSE115823", "GSE127813", "GSE53655", "GSE60424", "GSE64655", "sc_pbmc", "SDY67")
#blood_ref <- c("EPIC_BRef", "LM22_matrix_with_BC_CIBERSORTX.tsv", "Scaden_PBMC", "xCell2.0 - Blood")
#blood_tumor <- c("EPIC_BRef", "LM22_matrix_with_BC_CIBERSORTX.tsv", "Scaden_PBMC", "xCell2.0 - Blood")


xCell.data <- xCell::xCell.data

truths_dir <- "../xCell2.0/Kassandara_data/cell_values/"
mix_dir <- "../xCell2.0/Kassandara_data/expressions/"
ds <- gsub(".tsv", "", list.files(truths_dir))

score_xCell2 <- function(mix, sigs){
  mix_ranked <- singscore::rankGenes(mix)

  scores <- sapply(sigs, simplify = TRUE, function(sig){
    singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })

  if (is.list(scores)) {
    sigs <- sigs[-which(lengths(scores) == 0)]
    scores <- sapply(sigs, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
  }

  colnames(scores) <- names(sigs)
  rownames(scores) <- colnames(mix_ranked)
  scores <- as.data.frame(scores)
  scores <- scores %>%
    rownames_to_column(var = "sample") %>%
    pivot_longer(-sample, names_to = "signature", values_to = "score") %>%
    rowwise() %>%
    separate(signature, into = "celltype", sep = "#", remove = FALSE, extra = "drop") %>%
    group_by(sample, celltype) %>%
    summarise(score = median(score)) %>%
    pivot_wider(names_from = sample, values_from = score) %>%
    as_data_frame()
  xcell2.out <- data.frame(scores[,-1], row.names = scores$celltype, check.names = FALSE)

  return(xcell2.out)

}



getCorrelations <- function(dataset, celltypes2add = c("Monocytes", "Neutrophils", "NK cells", "cDC", "T-cells", "Tregs",
                                                       "Memory T-helpers", "CD4+ T-cells", "CD8+ T-cells", "B-cells",
                                                       "Fibroblasts", "Lymphocytes", "Cancer cells", "CD8+ T-cells PD1 low", "CD8+ T-cells PD1 high",
                                                       "Macrophages", "Macrophages M1", "Macrophages M2", "Endothelial cells", "Non plasma B-cells")){

  # Load mixture
  mix <- read.table(paste0("../xCell2.0/Kassandara_data/expressions/", dataset, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE)

  # Load truth
  truth <- read_tsv(paste0("../xCell2.0/Kassandara_data/cell_values/", dataset, ".tsv")) %>%
    dplyr::rename(celltype = 1) %>% # Fix for some datasets
    filter(!endsWith(celltype, "l)")) %>%
    filter(celltype != "Respiratory_cells") %>%
    filter(celltype != "Tumor KI67+") %>%
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
  xcell2_blood.out <- score_xCell2(mix, sigs = xcell2_blood_ref@filtered_signatures)
  xcell2_tumor.out <- score_xCell2(mix, sigs = xcell2_tumor_ref@filtered_signatures)
  xcell2_bp.out <- score_xCell2(mix, sigs = xcell2_bp_refsigs)


  # Run xCell
  xcell.out <- as.data.frame(xCell::xCellAnalysis(mix,rnaseq = FALSE))
  rownames(xcell.out) <- plyr::mapvalues(rownames(xcell.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

  # Load results from other methods
  kass.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/Kassandra/", dataset ,"_predicted_by_Kassandra.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(kass.out) <- plyr::mapvalues(rownames(kass.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  abis.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/ABIS/", dataset, "_predicted_by_ABIS.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(abis.out) <- plyr::mapvalues(rownames(abis.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  epic_bref.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/EPIC_BRef/", dataset, "_predicted_by_EPIC_BRef.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(epic_bref.out) <- plyr::mapvalues(rownames(epic_bref.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  epic_tref.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/EPIC_TRef/", dataset, "_predicted_by_EPIC_TRef.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(epic_tref.out) <- plyr::mapvalues(rownames(epic_tref.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  fardeep_abs.in <- paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/fardeep_absolute/", dataset, "_predicted_by_fardeep_absolute.tsv")
  if(file.exists(fardeep_abs.in)){
    fardeep_abs.out <- read.table(fardeep_abs.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(fardeep_abs.out) <- plyr::mapvalues(rownames(fardeep_abs.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    fardeep_abs.out <- NULL
  }
  fardeep_rel.in <- paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/fardeep_relative/", dataset, "_predicted_by_fardeep_relative.tsv")
  if(file.exists(fardeep_rel.in)){
    fardeep_rel.out <- read.table(fardeep_rel.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(fardeep_rel.out) <- plyr::mapvalues(rownames(fardeep_rel.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    fardeep_rel.out <- NULL
  }
  cbrx_hnsc.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/HNSC_matrix_with_BC_CIBERSORTX/", dataset, "_predicted_by_HNSC_matrix_with_BC_CIBERSORTX.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(cbrx_hnsc.out) <- plyr::mapvalues(rownames(cbrx_hnsc.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  cbrx_lm22.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/LM22_matrix_with_BC_CIBERSORTX/", dataset, "_predicted_by_LM22_matrix_with_BC_CIBERSORTX.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  cbrx_lm22.out <- cbrx_lm22.out[!rownames(cbrx_lm22.out) %in% c("Mast cells resting", "Mast cells activated", "NK cells resting", "NK cells activated",
                                                                 "Dendritic cells resting", "Dendritic cells activated"),] # Fix for some datasets
  rownames(cbrx_lm22.out) <- plyr::mapvalues(rownames(cbrx_lm22.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  quanti.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/quantiseq/", dataset, "_predicted_by_quantiseq.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(quanti.out) <- plyr::mapvalues(rownames(quanti.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  quanti_tumor.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/quantiseq_tumor/", dataset, "_predicted_by_quantiseq_tumor.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(quanti_tumor.out) <- plyr::mapvalues(rownames(quanti_tumor.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

  scaden.in <- paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/Scaden_PBMC/", dataset, "_predicted_by_Scaden_PBMC.tsv")
  if(file.exists(scaden.in)){
    scaden.out <- read.table(scaden.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(scaden.out) <- plyr::mapvalues(rownames(scaden.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    scaden.out <- NULL
  }


  all_methods <- c("xCell2.0 - Blood", "xCell2.0 - Tumor", "xCell2.0 - BlueprintEncode", "xCell", "Kassandara", "ABIS", "EPIC-BRef", "EPIC-TRef", "FARDEEP-Abs", "FARDEEP-Rel",
                   "CIBERSORTx-HNSC", "CIBERSORTx-LM22", "quanTIseq", "quanTIseq-T", "Scaden")

  all.out <- list(xcell2_blood.out, xcell2_tumor.out, xcell2_bp.out, xcell.out, kass.out, abis.out, epic_bref.out, epic_tref.out, fardeep_abs.out,
                  fardeep_rel.out, cbrx_hnsc.out, cbrx_lm22.out, quanti.out, quanti_tumor.out, scaden.out)
  names(all.out) <- all_methods

  all.out <- compact(all.out)
  shared_samples <- Reduce(intersect, lapply(all.out, names))
  shared_celltypes <- Reduce(intersect, lapply(all.out, rownames))
  shared_celltypes <- unique(c(shared_celltypes, celltypes2add))

  all.out <- lapply(all.out, function(x){x[shared_celltypes, shared_samples]})
  all.out <- lapply(all.out, function(x){cbind(celltype = shared_celltypes, x)})
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




cors.out.list <- vector("list", length(ds))
names(cors.out.list) <- ds

for (d in ds) {

#    if (d %in% ds[1:19]) {
#    next
#  }

  print(paste0(d, " --------------------------------------------------------- "))
  cors.out <- getCorrelations(dataset = d) %>%
    mutate(dataset = d)
  cors.out.list[[d]] <- cors.out
}


cors.out.list <- readRDS("/Users/almogang/Documents/xCell2_git/dev_scripts/cors.out.list_500genes_grubbs03_more_cts_grubbs05.rds")


all_cors.out <- bind_rows(cors.out.list)

all_cors.out <- all_cors.out %>%
  rowwise() %>%
  mutate(ds_type = ifelse(dataset %in% blood_ds, "Blood datasets", "Tumor datasets"))


ctoi <- "NK cells"

all_cors.out %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2.0 - Blood", "xCell2.0 - Tumor", "xCell2.0 - BlueprintEncode"), "yes", "no")) %>%
  ungroup() %>%
  filter(dataset != "Tonsils") %>%
  filter(cor_method == "Spearman") %>%
  filter(celltype == ctoi) %>%
  ggplot(., aes(x=method, y=cor_score, fill=is_xcell2)) +
  geom_boxplot(position = position_dodge(1), alpha = 0.6, outlier.shape = NA,  alpha = 0.5) +
  geom_point(aes(col=dataset), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  facet_grid(. ~ ds_type) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96")) +
  scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  labs(y = "Correlation coefficient (r)", title = ctoi) +
  theme_linedraw() +
  theme(strip.background =element_rect(fill="#8B2323"),
        plot.title = element_text(size=22),
        axis.ticks = element_line(linetype = "blank"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "#5E5E5E",
                                        linetype = "dashed"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 40, hjust=1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        strip.text.x = element_text(size = 16)) +
  labs(x = NULL, colour = NULL, fill = NULL)


cors.out.list[[1]] %>%
  filter(cor_method == "Spearman") %>%
  ggplot(., aes(x=method, y=cor_score)) +
  geom_boxplot(aes(fill=cor_method), position = position_dodge(1), alpha = 0.6, outlier.shape = NA) +
  geom_point(aes(col=celltype, group=cor_method), size=3,  position = position_jitterdodge(jitter.width = .3, dodge.width = 2)) +
  facet_grid(. ~ cor_method) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_fill_brewer(palette="Pastel1") +
  labs(y = "Correlation coefficient (r)") +
  theme(axis.ticks = element_line(linetype = "blank"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray95",
                                        linetype = "dashed"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 40, hjust=1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "gray97"),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        strip.text.x = element_text(size = 16)) +
  labs(x = NULL, colour = NULL, fill = NULL)






