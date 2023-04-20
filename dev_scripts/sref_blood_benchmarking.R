library(tidyverse)

sref.labels <- readRDS("/bigdata/almogangel/xCell2/data/sref_labels.rds")
sref.data <- readRDS("/bigdata/almogangel/xCell2/data/sref_data.rds")

# Load sref locations
sref.dic <- readxl::read_excel("/bigdata/almogangel/xCell2/dev_data/sref_cell_types_location.xlsx") %>%
  rowwise() %>%
  mutate(tissue.organ = str_split(tissue.organ, ", ")) %>%
  unnest(cols = c(tissue.organ))

# Use only blood samples
blood_samples <- sref.labels %>%
  left_join(sref.dic, by = c("ont", "label")) %>%
  filter(tissue.organ == "blood") %>%
  pull(sample)

sref.labels.blood <- sref.labels %>%
  filter(sample %in% blood_samples)

# Get blood validation datasets cell types
blood_ds <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813")

blood_celltypes <- c()
for (ds in blood_ds) {
  file.in <- read.table(paste0("/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/", ds, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  blood_celltypes <- unique(c(blood_celltypes, rownames(file.in)))
}

celltype_conversion_long <- read_tsv("dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

blood_cts <- plyr::mapvalues(blood_celltypes, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

# Missing cell types in reference
# TODO: Add a feature that find which cell types can form some of those cell types
missing_blood_cts <- blood_cts[!blood_cts %in% sref.labels.blood$label]


# Run xCell2.0
ref <- sref.data[,sref.labels.blood$sample]
labels <- sref.labels.blood

