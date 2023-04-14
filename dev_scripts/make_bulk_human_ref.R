library(tidyverse)

sref.df <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/sref_labels.rds")
sref.data <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/sref_data.rds")
kass.df <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/kass_labels.rds")
kass.data <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/kass_data.rds")

# TODO: !!! Remove shared samples from Mahmoud's reference (use Kassandra's)
# # Check for overlap between Mahmoud's and Kassandra's references
# sref_ds <- unique(unlist(strsplit(sref.df$dataset, ",")))
# kass_ds <- unique(kass.df$Dataset)
# shared_ds <- intersect(sref_ds, kass_ds)
#
# ds2remove <- !str_detect(sref.df$dataset, paste0("\\b", shared_ds, "\\b", collapse = "|"))
# sref.labels <- sref.df[ds2remove,]
# sref.data <- sref.df[,ds2remove]


human_bulk_labels.df <- rbind(sref.df, kass.df)

shared_genes <- intersect(rownames(sref.data), rownames(kass.data))
human_bulk_data.df <- cbind(sref.data[shared_genes,], kass.data[shared_genes,])

all(colnames(human_bulk_data.df) == human_bulk_labels.df$sample)

# Create S4 object for the new reference
setClass("xCell2 Bulk Human Reference", slots = list(
  labels = "data.frame",
  data = "matrix"))

xcell2_bulk_human_ref <- new("xCell2 Bulk Human Reference", labels = human_bulk_labels.df, data = human_bulk_data.df)

saveRDS(xcell2_bulk_human_ref, "/bigdata/almogangel/super_ref_for_xcell2/human_bulk_ref.rds")

