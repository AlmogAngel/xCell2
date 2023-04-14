# Make signatures for super reference

library(tidyverse)
source("/home/almogangel/xCell2_git/R/xCell2.R")

sref <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/human_bulk_ref.rds")

# # Use only main labels
# sref_labels_main <- sref@labels[,c(1,2,5,6)]
# colnames(sref_labels_main)[1:2] <- c("ont", "label")
# # xCell2GetLineage(sref_labels_main[,1:2], out_file = "/home/almogangel/xCell2_git/Data/sref_bulk_human_dependencies_main.tsv")

# Use main and fine labels
sref_labels_main <- sref@labels[,c(1,2,5,6)]
colnames(sref_labels_main)[1:2] <- c("ont", "label")

sref_labels_fine <- sref@labels[sref@labels$label.main != sref@labels$label.fine,c(3,4,5,6)]
colnames(sref_labels_fine)[1:2] <- c("ont", "label")

sref_labels <- rbind(sref_labels_main, sref_labels_fine)
# xCell2GetLineage(sref_labels[,1:2], out_file = "/home/almogangel/xCell2_git/Data/sref_bulk_human_dependencies.tsv")


ontology_file_checked <- "/home/almogangel/xCell2_git/Data/sref_bulk_human_dependencies_checked.tsv"

ref <- cbind(sref@data[,sref_labels_main$sample], sref@data[,sref_labels_fine$sample])
labels <- sref_labels
#dim(ref)
#dim(labels)
all(colnames(ref) == labels$sample)

sref_sigs <- xCell2Train(ref = ref, labels = labels, ontology_file_checked = ontology_file_checked, data_type = "rnaseq")

saveRDS(sref_sigs, "/home/almogangel/xCell2_git/Data/benchmarking_data/xcell2sigs_sref_main_fine.rds")
