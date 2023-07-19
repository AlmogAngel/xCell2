library(tidyverse)
library(Matrix)

setwd("/bigdata/almogangel/xCell2/")
source("R/xCell2Train.R")
source("R/xCell2GetLineage.R")

celltype_conversion_long <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))


######################## Single-cell RNA-seq references ---------------------------------------
############# Tabula Sapiens ---------------------------------------

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
  nest() %>%
  rowwise() %>%
  mutate(n_celltypes = length(unique(data$label)),
         celltypes = list(unique(data$label)))

# data = ts_labels.nested[3,]$data[[1]]
# tissue = ts_labels.nested[3,]$tissue[[1]]
makeRef <- function(data, tissue, ts_ref){
  print(tissue)

  xCell2GetLineage(data, out_file = paste0("/bigdata/almogangel/xCell2_data/dev_data/ts_", tissue, "_dependencies.tsv"))
  # View(read.table(paste0("/bigdata/almogangel/xCell2_data/dev_data/ts_", tissue, "_dependencies.tsv"), sep = "\t", header = TRUE))

  ref_tmp <- ts_ref[,data$sample]
  ref_tmp <- ref_tmp[Matrix::rowSums(ref_tmp) > 0,]

  ts_tissue_ref <- list(ref = ref_tmp,
                        labels = as.data.frame(data),
                        lineage_file = paste0("/bigdata/almogangel/xCell2_data/dev_data/ts_", tissue, "_dependencies.tsv"))

  ref_path <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/ts_", tissue, "_ref.rds")
  saveRDS(ts_tissue_ref, ref_path)

  return(ref_path)
}

# Generate references
ts_labels.nested %>%
  rowwise() %>%
  mutate(ref_path = makeRef(data, tissue, ts_ref))


###### Human (Tabula Sapiens) - blood -----------

# cytoVals - from run_validation.R
# val_celltypes <- unique(unname(unlist(sapply(cytoVals$truth$blood, rownames))))


# Subset blood related samples
ts_labels <- ts_labels[ts_labels$tissue %in% c("blood", "lymph node", "spleen", "thymus", "bone marrow", "inguinal lymph node"),]


# Use only shared cell types
shared <- intersect(val_celltypes, ts_labels$label)

# Use only cell type with at least 200 cells
shared <- names(which(sort(table(ts_labels[ts_labels$label %in% shared,]$label), decreasing = T) > 200))

# Subset cells
set.seed(123)

ts_labels_blood <- tibble(ts_labels) %>%
  filter(label %in% shared) %>%
  group_by(label, dataset) %>%
  sample_n(min(200, n())) %>%
  select(-tissue) %>%
  ungroup() %>%
  mutate(ont = as.character(ont),
         label = as.character(label),
         dataset = as.character(dataset))

labels <- data.frame(ts_labels_blood)
ref <- ts_ref[,labels$sample]
ref <- ref[Matrix::rowSums(ref) > 0,]


xCell2GetLineage(labels = labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv")

View(read.table("/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv", sep = "\t", header = TRUE))

ts_blood_ref <- list(ref = ref,
                       labels = labels,
                       lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv")
saveRDS(ts_blood_ref, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/ts_blood_ref.rds")



############# Pan-cancer ---------------------------------------

# Load counts matrix data
ref.mat.files <- list.files("/bigdata/almogangel/pan_cancer_ref/split_by_patient_matrices/", pattern = ".*RDS", full.names = T)
ref.mat.list <- lapply(ref.mat.files, function(path){
  ref.in <- readRDS(path)
  # Load only references with gene symbols:
  if ("CCDC7" %in% rownames(ref.in)) {
    ref.in
  }else{
    NA
  }
})

# Convert all matrices to dgCMatrix
ref.mat.list <- lapply(ref.mat.list, function(x){
  if ("matrix" %in% class(x)) {
    as(x, "dgCMatrix")
  }else{
    x
  }
})
ref.mat.list <- ref.mat.list[unlist(lapply(ref.mat.list, function(x){ncol(x) > 0}))]
ref.mat.list <- ref.mat.list[unlist(lapply(ref.mat.list, function(x){class(x) != "numeric"}))]


# Load metadata
meta <- readRDS("/bigdata/almogangel/pan_cancer_ref/harmonized_results_matrix_and_metadata_all_patients_all_cells.RDS")

# Define 10X datasets (using supplemantry table - 41467_2023_37353_MOESM11_ESM)
chromium10x_data <- c("Bi", "Chen", "Couterier", "Dong", "Gao", "Kim", "leader", "Lee", "Ma", "Peng", "Qian", "Slyper", "Wu", "Young")

meta.10x.confident <- meta %>%
  as_tibble() %>%
  # Remove non-10X data
  filter(dataset %in% chromium10x_data) %>%
  dplyr::select(cell_names, dataset, disease, primary_met, layer_1:layer_6, classification_confidence, scATOMIC_pred) %>%
  # Only cells with confident annotation
  filter(classification_confidence == "confident") %>%
  dplyr::select(-classification_confidence)

# Filter labels
labels2remove <- c("Tissue_Cell_Normal_or_Cancer", "unclassified_any_cell", "T_or_NK_lymphocyte", "macrophage_or_dendritic_cell", "unclassified_blood_cell",
                   "unclassified_normal_or_cancer_tissue", "NK or CD8 T cell", "Non GI Epithelial Cell", "unclassified_T_or_NK_cell", "unclassified_macrophage_or_DC",
                   "unclassified_B_cell_or_plasmablast", "GI Epithelial Cell", "Soft Tissue or Neuro Cancer Cell", "Ovarian/Endometrial/Kidney Cell", "unclassified_non_GI_epithelial_cell",
                   "unclassified_GI_epithelial_cell", "Breast/Lung/Prostate Cell", "Unclassified Soft Tissue or Neuro Cancer Cell", "Brain/Neuroblastoma Cancer Cell", "Any Cell",
                   "CD4 or CD8 T cell", "Non Stromal Cell", "Macrophage or Monocyte", "CD8 T or NK cell", "Ovarian/Endometrial/Kidney", "Normal Tissue Cell",
                   "B cell or Plasmablast", "Macrophage or Dendritic Cell", "Tfh/Th1 helper CD4+ T cells", "T or NK Cell", "Non Blood Cell",
                   "Biliary/Hepatic Cancer Cell", "Colorectal/Esophageal/Gastric Cell", "GI Tract Cell", "Breast/Lung/Prostate")

meta.10x.confident.filt <- meta.10x.confident %>%
  pivot_longer(c(layer_1:scATOMIC_pred), names_to = "layer", values_to = "label") %>%
  dplyr::select(-layer) %>%
  # Remove abritrary labels
  filter(!label %in% labels2remove) %>%
  unique() %>%
  mutate(label = plyr::mapvalues(label, from = celltype_conversion_long$all_labels, to = celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  group_by(label) %>%
  # Minimum 50 cells per label
  filter(n() >= 50) %>%
  ungroup()

# Subset and balance cells
set.seed(123)

meta.10x.confident.filt.sub <- meta.10x.confident.filt %>%
  group_by(label, dataset, disease, primary_met) %>%
  sample_n(min(1000, n())) %>%
  group_by(label, dataset, disease) %>%
  sample_n(min(1000, n())) %>%
  group_by(label, dataset) %>%
  sample_n(min(1000, n())) %>%
  group_by(label) %>%
  sample_n(min(1000, n()))

table(meta.10x.confident.filt.sub$label)
nrow(meta.10x.confident.filt.sub)
length(unique(meta.10x.confident.filt.sub$cell_names))

# Subset cells from reference data

ref.mat.list <- lapply(ref.mat.list, function(x){x[,colnames(x) %in% meta.10x.confident.filt.sub$cell_names]})
ref.mat.list <- ref.mat.list[unlist(lapply(ref.mat.list, function(x){class(x) != "numeric"}))]
ref.mat.list <- ref.mat.list[unlist(lapply(ref.mat.list, function(x){ncol(x) > 0}))]


# Combine references and make sure they share the same genes
shared_rownames <- Reduce(intersect, lapply(ref.mat.list, rownames))
ref.mat <- lapply(ref.mat.list, function(x) x[shared_rownames, ]) %>%
  do.call(cbind, .)

# Make sure every cell exist both in ref.mat and in the metadata
meta.10x.confident.filt.sub <- meta.10x.confident.filt.sub %>%
  filter(cell_names %in% colnames(ref.mat))
all(meta.10x.confident.filt.sub$cell_names %in% colnames(ref.mat))
table(meta.10x.confident.filt.sub$label)

ref.mat <- ref.mat[,meta.10x.confident.filt.sub$cell_names]
all(meta.10x.confident.filt.sub$cell_names == colnames(ref.mat))
meta.10x.confident.filt.sub$cell_names <- make.unique(meta.10x.confident.filt.sub$cell_names)
colnames(ref.mat) <- meta.10x.confident.filt.sub$cell_names

# Get dependencies
labels <- meta.10x.confident.filt.sub %>%
  mutate(ont = plyr::mapvalues(label, from = celltype_conversion_long$xCell2_labels, to = celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  dplyr::select(ont,  label, sample = cell_names, dataset)

xCell2GetLineage(labels = labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/sc_pan_cancer_dependencies.tsv")
lineage.file <- read.table("/bigdata/almogangel/xCell2_data/dev_data/sc_pan_cancer_dependencies.tsv", sep = "\t", header = TRUE)
lineage.file[lineage.file$label == "Exhausted CD8+ T cells",]$ancestors <- "CD8-positive, alpha-beta T cell"
lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             "Exhausted CD8+ T cells")
write.table(lineage.file, "/bigdata/almogangel/xCell2_data/dev_data/sc_pan_cancer_dependencies.tsv", sep = "\t", quote = F, row.names = F)



sc_pan_cancer_ref <- list(ref = ref.mat,
                          labels = as.data.frame(labels),
                          lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/sc_pan_cancer_dependencies.tsv")
saveRDS(sc_pan_cancer_ref, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/sc_pan_cancer_ref.rds")


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


lm22.ref <- read.table("/bigdata/almogangel/xCell2_data/dev_data/LM22-ref-sample.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
lm22.pheno <- read.table("/bigdata/almogangel/xCell2_data/dev_data/LM22-classes.txt", sep = "\t", header = FALSE, row.names = 1)
colnames(lm22.pheno) <- colnames(lm22.ref)
rownames(lm22.pheno) <- plyr::mapvalues(rownames(lm22.pheno), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels)


lm22.labels <- sapply(rownames(lm22.pheno), function(label){
  colnames(lm22.pheno)[which(lm22.pheno[label, ] == 1)]
}) %>%
  enframe(value = "sample", name = "label") %>%
  unnest(sample) %>%
  mutate(ont = NA, .before = everything()) %>%
  mutate(dataset = "LM22")

lm22.labels$ont <- plyr::mapvalues(lm22.labels$label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont)

lm22.ref <- lm22.ref[,lm22.labels$sample]

lm22.labels <- as.data.frame(lm22.labels)
xCell2GetLineage(lm22.labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")

lm22_ref <- list(ref = as.matrix(lm22.ref),
                       labels = lm22.labels,
                       lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")
saveRDS(lm22_ref, "/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")


# old
# lm22 <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22.rds")
# lm22.labels <- lm22$labels
#
# x <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/LM22.txt")
# x <- colnames(x)[-1]
# missing_cts <- unique(x[!x %in% lm22.labels])
#
#
# lm22.labels2 <- plyr::mapvalues(lm22.labels, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels)
# unique(lm22.labels2[lm22.labels2 %in% lm22.labels])
#
#
# lm22.labels <- lm22.labels2
# lm22.samples <- colnames(lm22$expr)
#
# lm22.dataset <- gsub(".CEL", "", lm22.samples)
# lm22.dataset <- gsub("_.*", "", lm22.dataset)
#
#
# x <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/lm22_datasets.txt")
#
# x$`Sample ID`[!x$`Sample ID` %in% lm22.dataset]
#
# unique(x$`Data set`[!x$`Data set` %in% lm22.dataset])
#
# which(startsWith(lm22.samples, "GSE4527"))
#
# lm22_labels <- tibble(ont = NA, label = lm22.labels, sample = lm22.samples, dataset = lm22.dataset)
# lm22_labels$ont <- plyr::mapvalues(lm22_labels$label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont)
# lm22_labels[!lm22_labels$ont %in% celltype_conversion_long$ont,]$ont <- NA
#
# all(lm22_labels$sample == colnames(lm22$expr))
#
#
# xCell2GetLineage(lm22_labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")
#
# lm22_ref <- list(ref = as.matrix(lm22$expr),
#                        labels = as.data.frame(lm22_labels),
#                        lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")
# saveRDS(lm22_ref, "/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")

