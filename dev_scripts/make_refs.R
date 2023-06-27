library(tidyverse)

setwd("/bigdata/almogangel/xCell2/")
source("R/xCell2Train.R")
source("R/xCell2GetLineage.R")

celltype_conversion_long <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))


######################## Single-cell RNA-seq references ---------------------------------------
############# Human (Tabula Sapiens)---------------------------------------

ts <- readRDS("/bigdata/almogangel/TabulaSapiens/tabula_sapiens.rds")
ts_ref <- ts@assays$RNA@counts
rownames(ts_ref) <- ts@assays$RNA@meta.features$feature_name

ts_labels <- ts@meta.data
ts_labels <- tibble(ts@meta.data) %>%
  mutate(sample = rownames(ts_labels)) %>%
  filter(assay == "10x 3' v3") %>%
  select(ont = cell_type_ontology_term_id, label = cell_type, sample, dataset = donor, tissue)

ts_ref <- ts_ref[,ts_labels$sample]


###### Human (Tabula Sapiens) - blood -----------

# cytoVals - from run_validation.R
val_celltypes <- unique(unname(unlist(sapply(cytoVals$truth$blood, rownames))))


# Subset blood related samples
ts_labels <- ts_labels[ts_labels$tissue %in% c("blood", "lymph node", "spleen", "thymus", "bone marrow", "inguinal lymph node"),]


# Use only shared cell types
shared <- intersect(val_celltypes, ts_labels$label)


# Subset cells
set.seed(123)

ts_labels_blood <- tibble(ts_labels) %>%
  filter(label %in% shared) %>%
  group_by(label, dataset) %>%
  sample_n(min(200, n())) %>%
  select(-tissue) %>%
  group_by(label) %>%
  filter(n() > 100) %>%

sort(table(ts_labels_blood$label), decreasing = T)

ts_labels_blood_main <- tibble(ts_labels) %>%
  filter(!sample %in% ts_labels_blood_fine$sample) %>%
  filter(label.main %in% shared_main) %>%
  group_by(label.main, dataset) %>%
  sample_n(min(50, n()))


labels <- data.frame(ont = c(ts_labels_blood_main$ont.main, ts_labels_blood_fine$ont.fine),
                        label = c(ts_labels_blood_main$label.main, ts_labels_blood_fine$label.fine),
                        sample = c(ts_labels_blood_main$sample, ts_labels_blood_fine$sample),
                        dataset = c(ts_labels_blood_main$dataset, ts_labels_blood_fine$dataset))

ref <- ts_ref[,ts_labels$sample]


xCell2GetLineage(labels = labels[,1:2], out_file = "/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv")

View(read.table("/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv", sep = "\t", header = TRUE))

ts_blood_ref <- list(ref = ref,
                       labels = as.data.frame(labels),
                       lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv")
saveRDS(ts_blood_ref, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/ts_blood_ref.rds")


########################  Bulk RNA-seq sorted cells references ---------------------------------------
############# Human ---------------------------------------
# Tumor ref (Kassandra) ----
# all_models_expr <- read_csv("../xCell2.0/Kassandara_data/all_models_expr.tsv")
all_models_expr <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_expr.tsv")
all_models_expr <- data.frame(all_models_expr[,-1], row.names = all_models_expr$Gene, check.names = F)
# all_models_annot <- read_csv("../xCell2.0/Kassandara_data/all_models_annot.tsv")
all_models_annot <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"
all(colnames(all_models_expr) == all_models_annot$Sample)


tumor_ref <- all_models_expr[,!is.na(all_models_annot$Tumor_model_annot)]
tumor_labels <- all_models_annot[!is.na(all_models_annot$Tumor_model_annot),]

tumor_labels <- tumor_labels %>%
  mutate(label = plyr::mapvalues(Tumor_model_annot, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  dplyr::select(ont, label, sample = Sample, dataset = Dataset)


ontology_file_checked <- "/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_dependencies_checked.tsv"

x <- read_tsv(ontology_file_checked)

unique(tumor_labels$label[which(!tumor_labels$label %in% x$label)])
tumor_labels[tumor_labels$label == "Non plasma B-cells",]$label <- "Non plasma B cell"
tumor_labels[tumor_labels$label == "CD8+ T-cells PD1 low",]$label <- "CD8+ T cell PD1 low"
tumor_labels[tumor_labels$label == "CD8+ T-cells PD1 high",]$label <- "CD8+ T cell PD1 high"
unique(x$label[which(!x$label %in% tumor_labels$label)])




kass_tumor_ref <- list(ref = as.matrix(tumor_ref),
     labels = as.data.frame(tumor_labels),
     lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_dependencies_checked.tsv")
kass_tumor_ref <- saveRDS(kass_tumor_ref, "/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_ref.rds")


ref <- as.matrix(tumor_ref)
labels <- as.data.frame(tumor_labels)
dim(ref)
dim(labels)

# Blood ref (Kassandra) ----
# all_models_expr <- read_csv("../xCell2.0/Kassandara_data/all_models_expr.tsv")
all_models_expr <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_expr.tsv")
all_models_expr <- data.frame(all_models_expr[,-1], row.names = all_models_expr$Gene, check.names = F)
# all_models_annot <- read_csv("../xCell2.0/Kassandara_data/all_models_annot.tsv")
all_models_annot <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"
all(colnames(all_models_expr) == all_models_annot$Sample)

blood_ref <- all_models_expr[,!is.na(all_models_annot$Blood_model_annot)]
blood_labels <- all_models_annot[!is.na(all_models_annot$Blood_model_annot),]
all(colnames(blood_ref) == blood_labels[!is.na(blood_labels$Blood_model_annot),]$Sample)

blood_labels <- blood_labels %>%
  mutate(label = plyr::mapvalues(Blood_model_annot, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  dplyr::select(ont, label, Sample, Dataset)

# Add lab data
laboratory_data_expressions <- read_tsv("/bigdata/almogangel/kassandra_data/sorted_cells/laboratory_data_expressions.tsv")
laboratory_data_expressions <- data.frame(laboratory_data_expressions[,-1], row.names = laboratory_data_expressions$Gene, check.names = F)
laboratory_data_annotation <- read_tsv("/bigdata/almogangel/kassandra_data/sorted_cells/laboratory_data_annotation.tsv")
colnames(laboratory_data_annotation)[1] <- "Sample"
all(colnames(laboratory_data_expressions) == laboratory_data_annotation$Sample)
laboratory_data_expressions <- laboratory_data_expressions[,laboratory_data_annotation$Sample]
all(colnames(laboratory_data_expressions) == laboratory_data_annotation$Sample)

lab_labels <- laboratory_data_annotation %>%
  mutate(label = plyr::mapvalues(Cell_type, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  mutate(Dataset = "Kassandra_lab") %>%
  dplyr::select(ont, label, Sample, Dataset)

blood_labels <- rbind(blood_labels, lab_labels)
blood_ref <- cbind(blood_ref, laboratory_data_expressions)
all(blood_labels$Sample == colnames(blood_ref))
colnames(blood_labels)[3:4] <- c("sample", "dataset")
# saveRDS(blood_labels, "Data/kass_blood_labels_with_lab.rds")

x <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/kass_blood_dependencies_checked.tsv")
View(x)



# ontology_file_checked <- "/xCell2.0/kass_blood_dependencies_checked.tsv"

kass_blood_ref <- list(ref = as.matrix(blood_ref),
               labels = as.data.frame(blood_labels),
               lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/kass_blood_dependencies_checked.tsv")
saveRDS(kass_blood_ref, "/bigdata/almogangel/xCell2_data/dev_data/kass_blood_ref.rds")



ref <- as.matrix(blood_ref)
labels <- as.data.frame(blood_labels)

# Subset data
sort(table(labels$label), decreasing = T)
celltypes2use <- c("CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "monocyte", "natural killer cell",
                   "neutrophil", "T cell", "B cell", "myeloid cell", "plasmacytoid dendritic cell, human", "granulocyte", "conventional dendritic cell", "plasma cell")
labels <- labels[labels$label %in% celltypes2use,]
ref <- ref[,labels$sample]

all(labels$sample == colnames(ref))



# BlueprintEncode ----

bp <- celldex::BlueprintEncodeData()

bp_ref <- as.matrix(bp@assays@data$logcounts)
bp_labels <- bp@colData

bp_labels <- bp_labels %>%
  as_tibble() %>%
  rename("label.ont" = "ont") %>%
  rowwise() %>%
  mutate(label = plyr::mapvalues(label.fine, from = celltype_conversion_long$all_labels, to = celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  ungroup() %>%
  mutate(sample = rownames(bp_labels),
         dataset = "Blueprint-Encode") %>%
  select(ont, label, sample, dataset)


x <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/bp_dependencies_checked.tsv")
x

all(bp_labels$sample == colnames(bp_ref))

bp_ref <- list(ref = as.matrix(bp_ref),
                       labels = as.data.frame(bp_labels),
                       lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/bp_dependencies_checked.tsv")
saveRDS(bp_ref, "/bigdata/almogangel/xCell2_data/dev_data/bp_ref.rds")




# Super reference -----

sref.labels <- readRDS("/bigdata/almogangel/xCell2_data/data/sref_labels.rds")
sref.data <- readRDS("/bigdata/almogangel/xCell2/data/sref_data.rds")

sref.metadata <- readxl::read_excel("/bigdata/almogangel/xCell2/data/sref_cell_types_metadata.xlsx") %>%
  rowwise() %>%
  mutate(tissue.organ = str_split(tissue.organ, ", ")) %>%
  unnest(cols = c(tissue.organ))

# Generate signatures for blood
blood_onts <- sref.metadata %>%
  filter(tissue.organ == "blood") %>%
  pull(ont) %>%
  unique()

sref.labels.blood <- sref.labels %>%
  filter(ont %in% blood_onts)

# Remove those cell types
sort(table(sref.labels.blood$label), decreasing = T)
celltype2remove <- c("cord blood hematopoietic stem cell", "peripheral blood mononuclear cell", "blood cell", "leukocyte")
sref.labels.blood <- sref.labels.blood %>%
  filter(!label %in% celltype2remove)

labels <- as.data.frame(sref.labels.blood)
ref <- sref.data[,labels$sample]
all(colnames(ref) == labels$sample)

saveRDS(labels, "/bigdata/almogangel/xCell2/dev_data/sref_blood_labels_bulk.rds")
saveRDS(ref, "/bigdata/almogangel/xCell2/dev_data/sref_blood_data_bulk.rds")



######################## Microarry sorted cells references ---------------------------------------
############# Human ---------------------------------------
# LM222 -----------


lm22 <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22.rds")
lm22.labels <- lm22$labels

x <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/LM22.txt")
x <- colnames(x)[-1]
missing_cts <- unique(x[!x %in% lm22.labels])


lm22.labels2 <- plyr::mapvalues(lm22.labels, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels)
unique(lm22.labels2[lm22.labels2 %in% lm22.labels])


lm22.labels <- lm22.labels2
lm22.samples <- colnames(lm22$expr)

lm22.dataset <- gsub(".CEL", "", lm22.samples)
lm22.dataset <- gsub("_.*", "", lm22.dataset)


x <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/lm22_datasets.txt")

x$`Sample ID`[!x$`Sample ID` %in% lm22.dataset]

unique(x$`Data set`[!x$`Data set` %in% lm22.dataset])

which(startsWith(lm22.samples, "GSE4527"))

lm22_labels <- tibble(ont = NA, label = lm22.labels, sample = lm22.samples, dataset = lm22.dataset)
lm22_labels$ont <- plyr::mapvalues(lm22_labels$label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont)
lm22_labels[!lm22_labels$ont %in% celltype_conversion_long$ont,]$ont <- NA

all(lm22_labels$sample == colnames(lm22$expr))


xCell2GetLineage(lm22_labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")

lm22_ref <- list(ref = as.matrix(lm22$expr),
                       labels = as.data.frame(lm22_labels),
                       lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")
saveRDS(lm22_ref, "/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")

