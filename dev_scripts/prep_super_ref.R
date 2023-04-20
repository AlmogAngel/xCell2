library(tidyverse)


# Load Super reference
sref <- readRDS("/bigdata/mahmoudy/humanCellTypeAtlas.rds")

# Make tibble
sref.labels <- tibble(ont = sref$cl.id, label = sref$lowest.label, sample = as.character(sref$samples), dataset = as.character(sref$series), info = sref$tissue_titles)
sref.labels[sref.labels$label == "B cell, CD19-positive",]$label <- "B cell" # Fix this label

# Add part of Kassandra's reference samples
kref <- readRDS("/bigdata/mahmoudy/kass_tpm.rds")
kref.labels <- readRDS("/bigdata/almogangel/xCell2/data/kass_labels_mahmoud.rds")





# # Load dictionary for main ont, main label and more
# sref.dic <- readxl::read_excel("/bigdata/almogangel/xCell2/data/sref_cell_types_dictionary_3groups.xlsx")
# # sref.dic %>%
# #   rowwise() %>%
# #   mutate(tissue.organ = str_split(tissue.organ, ";")) %>%
# #   unnest(cols = c(tissue.organ))
#
# sref.dic %>%
#   mutate(lowest_lable = label.fine) %>%
#   rowwise() %>%
#   mutate(lowest_lable = ifelse(is.na(lowest_lable), label.intermediate, lowest_lable)) %>%
#   mutate(lowest_lable = ifelse(is.na(lowest_lable), label.main, lowest_lable)) %>%
#   View()
#
# # Make tibble
# sref.tbl <- tibble(ont = sref$cl.id, label = sref$lowest.label, sample = as.character(sref$samples), dataset = as.character(sref$series), info = sref$tissue_titles)
#
# x <- sref.tbl %>%
#   rowwise() %>%
#   mutate(label_type = ifelse(label %in% sref.dic$label.main, "main",
#                              ifelse(label %in% sref.dic$ont.intermediate, "intermediate", "fine")))
#
#
# x <- sref.tbl %>%
#   rowwise() %>%
#   mutate(label.main = ifelse(label %in% sref.dic$label.main, sref.dic[sref.dic$label.main == label,]$label.main, NA)) %>%
#   mutate(label.intermediate = ifelse(label %in% sref.dic$label.intermediate, sref.dic[sref.dic$label.intermediate == label,]$label.intermediate, NA)) %>%
#   mutate(label.fine = ifelse(label %in% sref.dic$label.fine, sref.dic[sref.dic$label.fine == label,]$label.fine, NA)) %>%
#   mutate(label.intermediate = ifelse(label.fine %in% sref.dic$label.fine, sref.dic[sref.dic$label.fine == label.fine,]$label.intermediate, label.intermediate)) %>%
#   mutate(label.main = ifelse(label.fine %in% sref.dic$label.fine, sref.dic[sref.dic$label.fine == label.fine,]$label.main, label.main)) %>%
#   mutate(label.main = ifelse(label.intermediate %in% sref.dic$label.intermediate, sref.dic[sref.dic$label.intermediate == label.intermediate,]$label.main, label.main)) %>%
#   left_join(sref.dic, by = c(label = ))
#
#
#
# sref.tbl$ont.main <- plyr::mapvalues(sref.tbl$label.fine, sref.dic$label.fine, sref.dic$ont.main, warn_missing = FALSE)
# sref.tbl$label.main <- plyr::mapvalues(sref.tbl$label.fine, sref.dic$label.fine, sref.dic$label.main, warn_missing = FALSE)
# sref.tbl$family <- plyr::mapvalues(sref.tbl$label.fine, sref.dic$label.fine, sref.dic$family, warn_missing = FALSE)
# sref.tbl$tissue.organ <- plyr::mapvalues(sref.tbl$label.fine, sref.dic$label.fine, sref.dic$tissue.organ, warn_missing = FALSE)
#
# # Cell types to remove
# unique(sref.tbl[is.na(sref.tbl$ont.main),]$label.fine)
# sref.tbl <- sref.tbl[!is.na(sref.tbl$ont.main),]
#
# # Merge Super ref and Kassandra
# # Load Kassandra's reference
# # kref <- readRDS("/bigdata/mahmoudy/kass_tpm.rds")
# # ann <- readRDS("/bigdata/mahmoudy/kass_tpm_annotations.rds")
#
# # missing_kass_samples <- kann[!names(kann) %in% sref.tbl$sample]
#
# # TODO: Ask Mahmoud for the sample's ontology, dataset and info
# # tibble(ont.main = NA, label.main = NA, ont.fine = NA, label.fine = missing_kass_samples,
# #        sample = names(missing_kass_samples), dataset = NA, info = NA)
#
# # Save
# sref.labels <- as.data.frame(sref.tbl)
# sref.data <- as.matrix(sref$data[,sref.labels$sample])
# all(sref.labels$sample == colnames(sref.data))
#
# saveRDS(sref.labels, "/bigdata/almogangel/xCell2/data/sref_labels.rds")
# saveRDS(sref.data, "/bigdata/almogangel/xCell2/data/sref_data.rds")

