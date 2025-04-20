library(tidyverse)
library(Matrix)
library(limma)
library(ontologyIndex)

setwd("/bigdata/almogangel/xCell2/")
source("R/xCell2Train.R")
source("R/xCell2GetLineage.R")
cl <- ontoProc::getOnto(ontoname = "cellOnto", year_added = "2023")


getDependencies <- function(lineage_file_checked){
  ont <- read_tsv(lineage_file_checked, show_col_types = FALSE) %>%
    mutate_all(as.character)
  
  celltypes <- pull(ont[,2])
  celltypes <- gsub("_", "-", celltypes)
  dep_list <- vector(mode = "list", length = length(celltypes))
  names(dep_list) <- celltypes
  
  for (i in 1:nrow(ont)) {
    descendants <-  gsub("_", "-", strsplit(pull(ont[i,3]), ";")[[1]])
    descendants <- descendants[!is.na(descendants)]
    
    ancestors <-  gsub("_", "-", strsplit(pull(ont[i,4]), ";")[[1]])
    ancestors <- ancestors[!is.na(ancestors)]
    
    dep_list[[i]] <- list("descendants" = descendants, "ancestors" = ancestors)
    
  }
  
  return(dep_list)
}
get_ontology_id <- function(cell_type_name, ontology = cl) {
  # Find the ontology ID corresponding to the cell type name
  ontology_id <- names(ontology$name)[ontology$name == cell_type_name]
  
  if (length(ontology_id) > 0) {
    return(ontology_id)
  } else {
    return(NA) # Return NA if the cell type name is not found
  }
}
get_cell_type_name <- function(ontology_id, ontology = cl) {
  # Check if the ontology ID exists
  if (ontology_id %in% names(ontology$name)) {
    return(ontology$name[[ontology_id]])
  } else {
    return(NA) # Return NA if the ontology ID is not found
  }
}

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
         n_cells = nrow(data),
         celltypes = list(unique(data$label)))

# data = ts_labels.nested[3,]$data[[1]]
# tissue = ts_labels.nested[3,]$tissue[[1]]
makeRef <- function(data, tissue, ts_ref){
  print(tissue)
  
  xCell2GetLineage(data, out_file = paste0("/bigdata/almogangel/xCell2_data/references/human/ts/ts_", tissue, "_dependencies.tsv"))
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
cytoVals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
val_celltypes <- unique(unname(unlist(sapply(cytoVals$truth$blood, rownames))))


# Subset blood related samples
ts_labels <- ts_labels[ts_labels$tissue %in% c("Blood", "Lymph_Node", "Spleen", "Thymus", "Bone_Marrow", "Lymph_Node"),]


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

lineage.file <- read.table("/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv", sep = "\t", header = TRUE)
View(lineage.file)
lineage.file[lineage.file$label == "regulatory T cell",]$ancestors <- paste0(lineage.file[lineage.file$label == "regulatory T cell",]$ancestors, ";",
                                                                             "regulatory T cell;", "CD4-positive, alpha-beta T cell;", "CD8-positive, alpha-beta T cell")
lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             "regulatory T cell")
lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             "regulatory T cell")


# write_tsv(lineage.file, "/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv")


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
                   "Biliary/Hepatic Cancer Cell", "Colorectal/Esophageal/Gastric Cell", "GI Tract Cell", "Breast/Lung/Prostate", "Blood_Cell", "Billiary Cell", "MAIT cells", "Stromal Cell"
)

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

meta.10x.confident.filt.sub <- meta.10x.confident.filt.sub[meta.10x.confident.filt.sub$label != "blood cell",]


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
lineage.file[lineage.file$label == "Cancer Associated Fibroblasts",]$ancestors <- "fibroblast"
lineage.file[lineage.file$label == "fibroblast",]$descendants <- "Cancer Associated Fibroblasts"

lineage.file[lineage.file$label == "regulatory T cell",]$ancestors <- "CD4-positive, alpha-beta T cell;CD8-positive, alpha-beta T cell"
lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             "regulatory T cell")
lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             "regulatory T cell")

write.table(lineage.file, "/bigdata/almogangel/xCell2_data/dev_data/sc_pan_cancer_dependencies.tsv", sep = "\t", quote = F, row.names = F)



sc_pan_cancer_ref <- list(ref = ref.mat,
                          labels = as.data.frame(labels),
                          lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/sc_pan_cancer_dependencies.tsv")
saveRDS(sc_pan_cancer_ref, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/sc_pan_cancer_ref.rds")


###### MCA - blood -----------


load("/bigdata/almogangel/xCell2_data/mouse/MCDAA.rdata")
mca_labels <- read_csv("/bigdata/almogangel/xCell2_data/mouse/MCDAA_cellinfo.csv") %>%
  select(barcodes, cell_type)
colnames(mca_labels) <- c("sample", "label")

# Blood cells only
mca_blood_labels <- mca_labels[mca_labels$sample %in% rownames(data@meta.data[data@meta.data$tissue %in% c("PeripheralBlood", "Spleen", "BoneMarrow", "Thymus"),]),]
mca_blood_labels$dataset <- "MCA"


mca_blood_ref <- data[,mca_blood_labels$sample]
rm(data)
gc()


mca_blood_labels <- as_tibble(mca_blood_labels) %>%
  rowwise() %>%
  mutate(label = plyr::mapvalues(label, from = celltype_conversion_long$all_labels, to = celltype_conversion_long$xCell2_labels, warn_missing = FALSE))

table(mca_blood_labels$label)

mca_blood_labels <- mca_blood_labels %>%
  rowwise() %>%
  mutate(ont = get_ontology_id(label))


# Remove cells with low/high counts
n_counts <- colSums(mca_blood_ref)
n_genes <- colSums(mca_blood_ref@assays$RNA$counts > 0)

median(n_genes)
median(n_counts)

plot(sort(n_genes), sort(n_counts))
plot(sort(n_genes[which(n_genes <= 3000)]))
plot(sort(n_counts))

cells2use1 <- names(n_counts[n_counts >= 1500])
cells2use2 <- names(n_genes[n_genes >= 700])

cells2use <- intersect(cells2use1, cells2use2)

mca_blood_labels <- mca_blood_labels %>%
  filter(sample %in% cells2use)

ct2use <- sort(table(mca_blood_labels$label), decreasing = T)
ct2use <- names(ct2use[ct2use > 160])

mca_blood_labels <- mca_blood_labels %>%
  filter(label %in% ct2use)

table(mca_blood_labels$label)


# Subset cells
set.seed(123)

mca_blood_labels <- tibble(mca_blood_labels) %>%
  group_by(label) %>%
  sample_n(min(1000, n())) %>%
  select(ont, label, sample, dataset) %>%
  ungroup()

table(mca_blood_labels$label)




mca_blood_ref <- mca_blood_ref@assays$RNA$counts[,mca_blood_labels$sample]
all(colnames(mca_blood_ref) == mca_blood_labels$sample)


mca_blood_labels <- mca_blood_labels[!mca_blood_labels$label %in% c("Erythroblast", "Erythroid cell"),]
mca_blood_ref <- mca_blood_ref[,mca_blood_labels$sample]
all(colnames(mca_blood_ref) == mca_blood_labels$sample)


xCell2GetLineage(mca_blood_labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/mca_blood_dependencies.tsv")

lineage.file <- read.table("/bigdata/almogangel/xCell2_data/dev_data/mca_blood_dependencies.tsv", sep = "\t", header = TRUE)

lineage.file[lineage.file$label == "B cell",]$descendants <- "plasma cell"
lineage.file[lineage.file$label == "plasma cell",]$ancestors <- "B cell"

write.table(lineage.file, "/bigdata/almogangel/xCell2_data/dev_data/mca_blood_dependencies.tsv", sep = "\t", quote = F, row.names = F)




all(colnames(mca_blood_ref) == mca_blood_labels$sample)

mca_blood_ref <- mca_blood_ref[rowSums(mca_blood_ref) != 0,]


all(colnames(mca_blood_ref) == mca_blood_labels$sample)


mca_blood_ref <- list(ref = mca_blood_ref,
                      labels = as.data.frame(mca_blood_labels),
                      lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/mca_blood_dependencies.tsv")


saveRDS(mca_blood_ref, "/bigdata/almogangel/xCell2_data/dev_data/mca_blood_ref.rds")


# library(ExperimentHub)
# library(SingleCellExperiment)
#
#
# eh <- ExperimentHub()
# query(eh, "TabulaMurisData")
#
# tm.droplet <- eh[["EH1617"]]
# # tm.ss2 <- eh[["EH1618"]]
# # tm.droplet <- tm.ss2
#
# tm_ref <- tm.droplet@assays@data$counts
#
# tm_labels <- as_tibble(tm.droplet@colData) %>%
#   select(ont = cell_ontology_id, label = cell_ontology_class, sample = cell, dataset = mouse_id, tissue = tissue) %>%
#   drop_na()
#
# tm_ref <- tm_ref[,tm_labels$sample]
#
# all(colnames(tm_ref) == tm_labels$sample)
#
#
#
#
# # Subset blood related samples
# ts_labels <- ts_labels[ts_labels$tissue %in% c("blood", "lymph node", "Spleen", "Thymus", "Marrow", "inguinal lymph node"),]
#
#
# # Use only shared cell types
# shared <- intersect(val_celltypes, ts_labels$label)
#
# # Use only cell type with at least 200 cells
# shared <- names(which(sort(table(ts_labels[ts_labels$label %in% shared,]$label), decreasing = T) > 200))
#
# # Subset cells
# set.seed(123)
#
# ts_labels_blood <- tibble(ts_labels) %>%
#   filter(label %in% shared) %>%
#   group_by(label, dataset) %>%
#   sample_n(min(200, n())) %>%
#   select(-tissue) %>%
#   ungroup() %>%
#   mutate(ont = as.character(ont),
#          label = as.character(label),
#          dataset = as.character(dataset))
#
# labels <- data.frame(ts_labels_blood)
# ref <- ts_ref[,labels$sample]
# ref <- ref[Matrix::rowSums(ref) > 0,]
#
#
# xCell2GetLineage(labels = labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv")
#
# View(read.table("/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv", sep = "\t", header = TRUE))
#
# ts_blood_ref <- list(ref = ref,
#                      labels = labels,
#                      lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/ts_human_blood_dependencies.tsv")
# saveRDS(ts_blood_ref, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/ts_blood_ref.rds")





##### Tabula Muris blood ---------------
library(TabulaMurisData)


tm_blood <- TabulaMurisDroplet()

tm_blood_ref <- tm_blood@assays@data$counts
tissue2use <- tm_blood$tissue %in% c("Marrow", "Spleen", "Thymus")
ct2use <- tm_blood$cell_ontology_class %in% c("B cell", "T cell", "natural killer cell", "epithelial cell", "monocyte", "macrophage", "granulocyte")


cell2use <- tissue2use & ct2use
table(tm_blood$cell_ontology_class[cell2use])

tm_blood_ref <- tm_blood_ref[,cell2use]

tm_blood_labels <- tibble(ont = tm_blood$cell_ontology_id[cell2use], label = tm_blood$cell_ontology_class[cell2use],
                          sample = tm_blood$cell[cell2use], dataset = tm_blood$mouse_id[cell2use])

all(colnames(tm_blood_ref) == tm_blood_labels$sample)

xCell2GetLineage(tm_blood_labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/tm_blood_dependencies.tsv")

lineage.file <- read.table("/bigdata/almogangel/xCell2_data/dev_data/tm_blood_dependencies.tsv", sep = "\t", header = TRUE)

tm_blood_ref <- list(ref = tm_blood_ref,
                     labels = as.data.frame(tm_blood_labels),
                     lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/tm_blood_dependencies.tsv")


saveRDS(tm_blood_ref, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/tm_blood_ref.rds")



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
lineage.file <- read_tsv(ontology_file_checked)

View(lineage.file)
lineage.file[lineage.file$label == "regulatory T cell",]$ancestors <- paste0(lineage.file[lineage.file$label == "regulatory T cell",]$ancestors, ";",
                                                                             "CD8-positive, alpha-beta T cell")
lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants <- paste0("regulatory T cell")
lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             "regulatory T cell")

lineage.file[lineage.file$label == "helper T cell",]$ancestors <- paste0(lineage.file[lineage.file$label == "helper T cell",]$ancestors, ";",
                                                                         "CD4-positive, alpha-beta T cell")
lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             lineage.file[lineage.file$label == "helper T cell",]$descendants, ";",
                                                                                             "helper T cell")

lineage.file[lineage.file$label == "T follicular helper cell",]$ancestors <- paste0(lineage.file[lineage.file$label == "T follicular helper cell",]$ancestors, ";",
                                                                                    "CD4-positive, alpha-beta T cell")

lineage.file[lineage.file$label == "T-helper 1 cell",]$ancestors <- paste0(lineage.file[lineage.file$label == "T-helper 1 cell",]$ancestors, ";",
                                                                           "CD4-positive, alpha-beta T cell")
lineage.file[lineage.file$label == "T-helper 2 cell",]$ancestors <- paste0(lineage.file[lineage.file$label == "T-helper 2 cell",]$ancestors, ";",
                                                                           "CD4-positive, alpha-beta T cell")
lineage.file[lineage.file$label == "T-helper 17 cell",]$ancestors <- paste0(lineage.file[lineage.file$label == "T-helper 17 cell",]$ancestors, ";",
                                                                            "CD4-positive, alpha-beta T cell")


# write_tsv(lineage.file, "/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_dependencies_checked.tsv")



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




# rm_ct <- blood_labels$label != "myeloid cell"
# blood_labels <- blood_labels[rm_ct,]
# blood_ref <- blood_ref[,blood_labels$sample]
# all(blood_labels$sample == colnames(blood_ref))

lineage.file <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/kass_blood_dependencies_checked.tsv")
View(lineage.file)

# lineage.file[lineage.file$label == "T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "T cell",]$descendants, ";",
#                                                                                              "CD8-positive, alpha-beta memory T cell")
# lineage.file[lineage.file$label == "natural killer cell",]$descendants <- paste0(lineage.file[lineage.file$label == "natural killer cell",]$descendants, ";",
#                                                                     "Regulatory NK cells")
# lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants, ";",
#                                                                                  "CD8-positive, alpha-beta memory T cell")
# lineage.file[lineage.file$label == "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",]$ancestors <- paste0(lineage.file[lineage.file$label == "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",]$ancestors, ";",
#                                                                                              "effector memory CD8-positive, alpha-beta T cell")




write.table(lineage.file, "/bigdata/almogangel/xCell2_data/dev_data/kass_blood_dependencies_checked.tsv", sep = "\t", quote = F, row.names = F)



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

# Anti log2 data
bp_ref <- 2^bp_ref
bp_ref <- bp_ref-1

bp_labels <- bp@colData

bp_labels <- bp_labels %>%
  as_tibble() %>%
  rename("ont" = "label.ont") %>%
  rowwise() %>%
  mutate(label = plyr::mapvalues(label.fine, from = celltype_conversion_long$all_labels, to = celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  ungroup() %>%
  mutate(sample = rownames(bp_labels),
         dataset = "Blueprint-Encode") %>%
  select(ont, label, sample, dataset)




lineage.file <- read.table("/bigdata/almogangel/xCell2_data/dev_data/bp_dependencies_checked.tsv", sep = "\t", header = TRUE)
View(lineage.file)
lineage.file[lineage.file$label == "regulatory T cell",]$ancestors <- paste0("CD4-positive, alpha-beta T cell;", "CD8-positive, alpha-beta T cell")
lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD4-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             "regulatory T cell")
lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants <- paste0(lineage.file[lineage.file$label == "CD8-positive, alpha-beta T cell",]$descendants, ";",
                                                                                             "regulatory T cell")


# write_tsv(lineage.file, "/bigdata/almogangel/xCell2_data/dev_data/bp_dependencies_checked.tsv")


all(bp_labels$sample == colnames(bp_ref))

bp_ref <- list(ref = as.matrix(bp_ref),
               labels = as.data.frame(bp_labels),
               lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/bp_dependencies_checked.tsv")
saveRDS(bp_ref, "/bigdata/almogangel/xCell2_data/dev_data/bp_ref.rds")




# DICE ----

dice <- celldex::DatabaseImmuneCellExpressionData()

dice_ref <- as.matrix(dice@assays@data$logcounts)
dice_labels <- as.data.frame(dice@colData)
dice_labels$sample <- make.unique(rownames(dice_labels))
colnames(dice_ref) <- dice_labels$sample

cl <- ontoProc::getOnto(ontoname = "cellOnto", year_added = "2023")

dice_labels$label.fine <- unname(cl$name[dice_labels$label.ont])

dice_labels <- dice_labels %>%
  as_tibble() %>%
  dplyr::rename("ont" = "label.ont",
                "label" = "label.fine") %>%
  ungroup() %>%
  mutate(dataset = "DICE") %>%
  select(ont, label, sample, dataset)

# Switch labels
dice_labels[dice_labels$label == "activated CD4-positive, alpha-beta T cell", ]$ont <- "CL:0000624"
dice_labels[dice_labels$label == "activated CD4-positive, alpha-beta T cell", ]$label <- "CD4-positive, alpha-beta T cell"
dice_labels[dice_labels$label == "activated CD8-positive, alpha-beta T cell", ]$ont <- "CL:0000625"
dice_labels[dice_labels$label == "activated CD8-positive, alpha-beta T cell", ]$label <- "CD8-positive, alpha-beta T cell"

# Remove labels
labels2remove <- !dice_labels$label %in% c("CD14-positive, CD16-negative classical monocyte", "CD14-low, CD16-positive monocyte",
                                           "naive CCR4-positive regulatory T cell")
dice_labels <- dice_labels[labels2remove,]
dice_ref <- dice_ref[,labels2remove]

xCell2GetLineage(labels = dice_labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/dice_dependencies.tsv")


all(dice_labels$sample == colnames(dice_ref))

dice_ref <- list(ref = as.matrix(dice_ref),
                 labels = as.data.frame(dice_labels),
                 lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/dice_dependencies.tsv")

saveRDS(dice_ref, "/bigdata/almogangel/xCell2_data/dev_data/dice_ref.rds")





############# Mouse ---------------------------------------
# ImmGenData ----

igd <- celldex::ImmGenData()

igd_ref <- as.matrix(igd@assays@data$logcounts)

# Anti log2 data
igd_ref <- 2^igd_ref
igd_ref <- igd_ref-1


igd_ref <- limma::normalizeBetweenArrays(igd_ref)
igd_ref <- log2(igd_ref + 1)


igd_labels <- igd@colData


# igd_labels <- as_tibble(igd_labels) %>%
#   rowwise() %>%
#   mutate(ont.label = get_cell_type_name(label.ont))


igd_labels <- igd_labels %>%
  as_tibble() %>%
  dplyr::rename("ont" = "label.ont") %>%
  rowwise() %>%
  mutate(label = plyr::mapvalues(label.main, from = celltype_conversion_long$all_labels, to = celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  ungroup() %>%
  mutate(sample = rownames(igd_labels),
         dataset = "ImmGenData") %>%
  select(ont, label, sample, dataset)

all(colnames(igd_ref) == igd_labels$sample)


igd_labels_fine <- igd@colData
igd_labels_fine <- igd_labels_fine %>%
  as_tibble() %>%
  dplyr::rename("ont" = "label.ont") %>%
  rowwise() %>%
  mutate(label = plyr::mapvalues(label.fine, from = celltype_conversion_long$all_labels, to = celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  ungroup() %>%
  mutate(sample = rownames(igd_labels_fine),
         dataset = "ImmGenData") %>%
  select(ont, label, sample, dataset)
igd_labels_fine[grepl("CD4", igd_labels_fine$label),]$label <- "CD4-positive, alpha-beta T cell"
igd_labels_fine[grepl("CD8", igd_labels_fine$label),]$label <- "CD8-positive, alpha-beta T cell"

samples2dup <- igd_labels_fine[igd_labels_fine$label %in% c("CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"),]$sample
labels_dup <- igd_labels_fine[igd_labels_fine$sample %in% samples2dup,]



igd_ref_dup <- igd_ref[,labels_dup$sample]
labels_dup$sample <- paste0(labels_dup$sample, "_dup")

colnames(igd_ref_dup) <- labels_dup$sample

igd_ref <- cbind(igd_ref, igd_ref_dup)
igd_labels <- rbind(igd_labels, labels_dup)

all(colnames(igd_ref) == igd_labels$sample)

igd_labels <- igd_labels %>%
  rowwise() %>%
  mutate(ont = get_ontology_id(label))

igd_labels <- igd_labels %>%
  drop_na()

igd_ref <- igd_ref[,igd_labels$sample]
all(colnames(igd_ref) == igd_labels$sample)


igd_labels <- igd_labels[igd_labels$label != "pro-B cell",]
igd_ref <- igd_ref[,igd_labels$sample]
all(colnames(igd_ref) == igd_labels$sample)

igd_labels <- as.data.frame(igd_labels)

xCell2GetLineage(igd_labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/igd_dependencies.tsv")
dep_list <- getDependencies("/bigdata/almogangel/xCell2_data/dev_data/igd_dependencies.tsv")

x <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/igd_dependencies.tsv")
x


igd_ref <- list(ref = as.matrix(igd_ref),
                labels = as.data.frame(igd_labels),
                lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/igd_dependencies.tsv")
saveRDS(igd_ref, "/bigdata/almogangel/xCell2_data/dev_data/igd_ref.rds")





# MouseRNAseqData ----

mouse_rnaseq_data <- celldex::MouseRNAseqData()

mouse_rnaseq_data_ref <- as.matrix(mouse_rnaseq_data@assays@data$logcounts)

# Anti log2 data
mouse_rnaseq_data_ref <- 2^mouse_rnaseq_data_ref
mouse_rnaseq_data_ref <- mouse_rnaseq_data_ref-1

mouse_rnaseq_data_labels <- mouse_rnaseq_data@colData


# igd_labels <- as_tibble(igd_labels) %>%
#   rowwise() %>%
#   mutate(ont.label = get_cell_type_name(label.ont))


mouse_rnaseq_data_labels <- mouse_rnaseq_data_labels %>%
  as_tibble() %>%
  dplyr::rename("ont" = "label.ont") %>%
  rowwise() %>%
  mutate(label = plyr::mapvalues(label.main, from = celltype_conversion_long$all_labels, to = celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  ungroup() %>%
  mutate(sample = rownames(mouse_rnaseq_data_labels),
         dataset = "MouseRNAseqData") %>%
  select(ont, label, sample, dataset)

all(colnames(mouse_rnaseq_data_ref) == mouse_rnaseq_data_labels$sample)



xCell2GetLineage(mouse_rnaseq_data_labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/mouse_rnaseq_data_dependencies.tsv")
dep_list <- getDependencies("/bigdata/almogangel/xCell2_data/dev_data/mouse_rnaseq_data_dependencies.tsv")

x <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/mouse_rnaseq_data_dependencies.tsv")
x


mouse_rnaseq_data_ref <- list(ref = as.matrix(mouse_rnaseq_data_ref),
                              labels = as.data.frame(mouse_rnaseq_data_labels),
                              lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/mouse_rnaseq_data_dependencies.tsv")
saveRDS(mouse_rnaseq_data_ref, "/bigdata/almogangel/xCell2_data/benchmarking_data/references/mouse_rnaseq_data_ref.rds")






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

# Remove cell types
# unique(lm22.labels$label)
# lm22.labels[lm22.labels$label == "T cells CD4 memory activated",]$label <- "CD4-positive, alpha-beta memory T cell"
# lm22.labels[lm22.labels$label == "T cells CD4 memory resting",]$label <- "CD4-positive, alpha-beta memory T cell"
# lm22.labels[lm22.labels$label == "Mast cells resting",]$label <- "mast cell"
# lm22.labels[lm22.labels$label == "Mast cells activated",]$label <- "mast cell"
# lm22.labels[lm22.labels$label == "NK cells activated",]$label <- "natural killer cell"
# lm22.labels[lm22.labels$label == "NK cells resting",]$label <- "natural killer cell"
# lm22.labels[lm22.labels$label == "Dendritic cells activated",]$label <- "dendritic cell"
# lm22.labels[lm22.labels$label == "Dendritic cells resting",]$label <- "dendritic cell"
# lm22.labels[lm22.labels$label == "inflammatory macrophage",]$label <- "macrophage"
# lm22.labels[lm22.labels$label == "alternatively activated macrophage",]$label <- "macrophage"
# lm22.labels$ont <- plyr::mapvalues(lm22.labels$label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont)

# samples2remove <- !lm22.labels$label %in% c("T cells CD4 memory activated", "T cells CD4 memory resting",
#                                            "Mast cells resting", "Mast cells activated",
#                                            "NK cells activated", "NK cells resting",
#                                            "Dendritic cells activated", "Dendritic cells resting")
# lm22.labels <- lm22.labels[samples2remove,]
# lm22.ref <- lm22.ref[,samples2remove]
# all(lm22.labels$sample == colnames(lm22.ref))

xCell2GetLineage(lm22.labels, out_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")

lineage.file <- read.table("/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv", sep = "\t", header = TRUE)




# lm22_ref <- list(ref = as.matrix(lm22.ref),
#                        labels = lm22.labels,
#                        lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")
#
# saveRDS(lm22_ref, "/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")


lm22.ref.rma <- backgroundCorrect(lm22.ref, method='normexp')
lm22.ref.rma <- normalizeBetweenArrays(lm22.ref.rma)
lm22.ref.rma <- log2(lm22.ref.rma + 1)

lm22_ref_rma <- list(ref = as.matrix(lm22.ref.rma),
                     labels = lm22.labels,
                     lineage_file = "/bigdata/almogangel/xCell2_data/dev_data/lm22_dependencies.tsv")

saveRDS(lm22_ref_rma, "/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")

