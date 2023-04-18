library(tidyverse)
library(ontologyIndex)

cl <- ontoProc::getCellOnto()
celltype_conversion_long <- read_tsv("dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))


# Load Super reference
sref <- readRDS("/bigdata/mahmoudy/humanCellTypeAtlas.rds")

# Load Kassandra's reference
kref <- readRDS("/bigdata/mahmoudy/kass_tpm.rds")
kann <- readRDS("/bigdata/mahmoudy/kass_tpm_annotations.rds")

# Make tibble
sref.tbl <- tibble(ont.main = NA, label.main = NA, ont.fine = sref$cl.id, label.fine = sref$lowest.label,
       sample = as.character(sref$samples), dataset = as.character(sref$series), info = sref$tissue_titles,
       family = NA, tissue.organ = NA)

# Load dictionary for main ont, main label and more
sref.dic <- read_tsv("/bigdata/almogangel/xCell2/data/sref_cell_types_dictionary.txt")
sref.tbl$ont.main <- plyr::mapvalues(sref.tbl$label.fine, sref.dic$label.fine, sref.dic$ont.main, warn_missing = FALSE)
sref.tbl$label.main <- plyr::mapvalues(sref.tbl$label.fine, sref.dic$label.fine, sref.dic$label.main, warn_missing = FALSE)
sref.tbl$family <- plyr::mapvalues(sref.tbl$label.fine, sref.dic$label.fine, sref.dic$family, warn_missing = FALSE)
sref.tbl$tissue.organ <- plyr::mapvalues(sref.tbl$label.fine, sref.dic$label.fine, sref.dic$tissue.organ, warn_missing = FALSE)

# Cell types to remove
unique(sref.tbl[is.na(sref.tbl$ont.main),]$label.fine)
sref.tbl <- sref.tbl[!is.na(sref.tbl$ont.main),]

# Merge Super ref and Kassandra
# missing_kass_samples <- kann[!names(kann) %in% sref.tbl$sample]

# TODO: Ask Mahmoud for the sample's ontology, dataset and info
# tibble(ont.main = NA, label.main = NA, ont.fine = NA, label.fine = missing_kass_samples,
#        sample = names(missing_kass_samples), dataset = NA, info = NA)

# Save
sref.labels <- as.data.frame(sref.tbl)
sref.data <- as.matrix(sref$data[,sref.labels$sample])
all(sref.labels$sample == colnames(sref.data))

saveRDS(sref.labels, "/bigdata/almogangel/xCell2/data/sref_labels.rds")
saveRDS(sref.data, "/bigdata/almogangel/xCell2/data/sref_data.rds")

