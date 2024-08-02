library(tidyverse)
library(xCell2)


celltype_conversion_long <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

min_ct_cells_value <- 50


# Human ---------------------------------------
# Tabula Sapiens ----

# Load TS
ts <- readRDS("/bigdata/almogangel/TabulaSapiens/tabula_sapiens.rds")
ts_ref <- ts@assays$RNA@counts
rownames(ts_ref) <- ts@assays$RNA@meta.features$feature_name

# Only 10X data
ts_labels <- as_tibble(ts@meta.data, rownames = "sample") %>%
  filter(assay == "10x 3' v3") %>%
  select(ont = cell_type_ontology_term_id, label = cell_type, sample, dataset = donor, tissue = tissue_in_publication)

ts_ref <- ts_ref[,ts_labels$sample]
rm(ts)
gc()

# Nest by organ/tissue
ts_labels.nested <- tibble(ts_labels) %>%
  mutate_if(is.factor, as.character) %>%
  group_by(tissue) %>%
  nest(.key = "labels") %>%
  rowwise() %>%
  mutate(n_celltypes = length(unique(labels$label)),
         n_cells = nrow(labels)) %>%
  mutate(n_celltypes_passed = sum(table(labels$label) >= min_ct_cells_value),
         n_cells_passed = sum(table(labels$label)[which(table(labels$label) >= min_ct_cells_value)]))

# labels = ts_labels.nested[3,]$labels[[1]]
# tissue = ts_labels.nested[3,]$tissue[[1]]
# ts_ref = ts_labels.nested[3,]$tissue[[1]]

makeRef <- function(labels, tissue, TSref = ts_ref){
  print(tissue)

  labels <- as.data.frame(labels)
  ref <- TSref[,labels$sample]
  ref <- ref[Matrix::rowSums(ref) > 0,]

  ct2use <- names(which(table(labels$label) >= min_ct_cells_value)) # Must have at least min_ct_cells_value cells
  labels <- labels[labels$label %in%ct2use, ]
  ref <- ref[, labels$sample]


  ref <- xCell2::xCell2Train(ref = as.matrix(ref), labels = as.data.frame(labels), ref_type = "sc",
                             num_threads = 40, use_ontology = TRUE)

  saveRDS(ref, paste0("/bigdata/almogangel/xCell2_data/references/human/ts/ts_", tissue, ".xcell2object.rds"))

  print("Done.")
}


for (i in 1:nrow(ts_labels.nested)) {

  file <- paste0("/bigdata/almogangel/xCell2_data/references/human/ts/ts_", ts_labels.nested$tissue[i], ".xcell2object.rds")
  if (file.exists(file)) {
    print("Next...")
    next
  }
  makeRef(labels = ts_labels.nested$labels[[i]], tissue = ts_labels.nested$tissue[i])
  gc()

}



# Tabula Muris ----

library(ExperimentHub)
library(SingleCellExperiment)


eh <- ExperimentHub()
query(eh, "TabulaMurisData")

tm.droplet <- eh[["EH1617"]]
# tm.ss2 <- eh[["EH1618"]]
# tm.droplet <- tm.ss2


tm_ref <- tm.droplet@assays@data$counts


tm_labels <- as_tibble(tm.droplet@colData) %>%
  select(ont = cell_ontology_id, label = cell_ontology_class, sample = cell, dataset = mouse_id, tissue = tissue) %>%
  drop_na()

tm_ref <- tm_ref[,tm_labels$sample]

all(colnames(tm_ref) == tm_labels$sample)

labels <- tm_labels

tm_labels.nested <- tibble(tm_labels) %>%
  mutate_if(is.factor, as.character) %>%
  group_by(tissue) %>%
  nest(.key = "labels") %>%
  rowwise() %>%
  mutate(n_celltypes = length(unique(labels$label)),
         n_cells = nrow(labels)) %>%
  mutate(n_celltypes_passed = sum(table(labels$label) >= min_ct_cells_value),
         n_cells_passed = sum(table(labels$label)[which(table(labels$label) >= min_ct_cells_value)]))


