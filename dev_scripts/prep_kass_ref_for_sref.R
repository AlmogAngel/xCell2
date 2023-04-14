library(tidyverse)
library(ontologyIndex)

cl <- ontoProc::getCellOnto()
celltype_conversion_long <- read_tsv("Data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# Load Kassandra's data
all_models_expr <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_expr.tsv")
all_models_expr <- data.frame(all_models_expr[,-1], row.names = all_models_expr$Gene, check.names = F)
all_models_annot <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"
laboratory_data_expressions <- read_tsv("/bigdata/almogangel/kassandra_data/sorted_cells/laboratory_data_expressions.tsv")
laboratory_data_expressions <- data.frame(laboratory_data_expressions[,-1], row.names = laboratory_data_expressions$Gene, check.names = F)
laboratory_data_annotation <- read_tsv("/bigdata/almogangel/kassandra_data/sorted_cells/laboratory_data_annotation.tsv")
colnames(laboratory_data_annotation)[1] <- "Sample"
laboratory_data_expressions <- laboratory_data_expressions[,laboratory_data_annotation$Sample]


# Blood samples first
all_models_annot_blood <- all_models_annot[!is.na(all_models_annot$Blood_model_annot),]
kass_blood.tbl <- tibble(ont.main = NA, label.main = NA, ont.fine = NA, label.fine = all_models_annot_blood$Blood_model_annot,
                   sample = all_models_annot_blood$Sample, dataset = all_models_annot_blood$Dataset, info = "Kassanda blood model")

kass_blood_lab.tbl <- tibble(ont.main = NA, label.main = NA, ont.fine = NA, label.fine = laboratory_data_annotation$Cell_type,
                         sample = laboratory_data_annotation$Sample, dataset = "Kassandra", info = "Kassanda blood lab")

kass_blood.tbl <- rbind(kass_blood.tbl, kass_blood_lab.tbl)

# Tumor samples
all_models_annot_tumor <- all_models_annot[!all_models_annot$Sample %in% kass_blood.tbl$sample,]
kass_tumor.tbl <- tibble(ont.main = NA, label.main = NA, ont.fine = NA, label.fine = all_models_annot_tumor$Tumor_model_annot,
                         sample = all_models_annot_tumor$Sample, dataset = all_models_annot_tumor$Dataset, info = "Kassanda tumor model")

# Merge all samples
kass_ref.tbl <- rbind(kass_blood.tbl, kass_tumor.tbl)

# Add labels and ontology
kass_ref.tbl <- kass_ref.tbl %>%
  mutate(label.fine = plyr::mapvalues(label.fine, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont.fine = plyr::mapvalues(label.fine, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  rowwise() %>%
  mutate(label.fine = ifelse(is.na(unname(cl$name[ont.fine])), label.fine, unname(cl$name[ont.fine])))


# Manually add main ontology and label
kass_ref.tbl %>%
  select(ont.fine, label.fine) %>%
  unique() %>%
  print(n=Inf)

kass_ref.tbl[kass_ref.tbl$label.fine == "CD4-positive, alpha-beta T cell",]$ont.main <- "CL:0000624"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD4-positive, alpha-beta T cell",]$label.main <- "CD4-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8-positive, alpha-beta T cell",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8-positive, alpha-beta T cell",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "granulocyte",]$ont.main <- "CL:0000094"
kass_ref.tbl[kass_ref.tbl$label.fine == "granulocyte",]$label.main <- "granulocyte"
kass_ref.tbl[kass_ref.tbl$label.fine == "natural killer cell",]$ont.main <- "CL:0000623"
kass_ref.tbl[kass_ref.tbl$label.fine == "natural killer cell",]$label.main <- "natural killer cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "monocyte",]$ont.main <- "CL:0000576"
kass_ref.tbl[kass_ref.tbl$label.fine == "monocyte",]$label.main <- "monocyte"
kass_ref.tbl[kass_ref.tbl$label.fine == "neutrophil",]$ont.main <- "CL:0000775"
kass_ref.tbl[kass_ref.tbl$label.fine == "neutrophil",]$label.main <- "neutrophil"
kass_ref.tbl[kass_ref.tbl$label.fine == "T cell",]$ont.main <- "CL:0000084"
kass_ref.tbl[kass_ref.tbl$label.fine == "T cell",]$label.main <- "T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "Non plasma B-cells",]$ont.main <- "CL:0000236"
kass_ref.tbl[kass_ref.tbl$label.fine == "Non plasma B-cells",]$label.main <- "B cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "B cell",]$ont.main <- "CL:0000236"
kass_ref.tbl[kass_ref.tbl$label.fine == "B cell",]$label.main <- "B cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "plasmacytoid dendritic cell, human",]$ont.main <- "CL:0000451"
kass_ref.tbl[kass_ref.tbl$label.fine == "plasmacytoid dendritic cell, human",]$label.main <- "dendritic cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "myeloid cell",]$ont.main <- "CL:0000763"
kass_ref.tbl[kass_ref.tbl$label.fine == "myeloid cell",]$label.main <- "myeloid cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "unswitched memory B cell",]$ont.main <- "CL:0000236"
kass_ref.tbl[kass_ref.tbl$label.fine == "unswitched memory B cell",]$label.main <- "B cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "non-classical monocyte",]$ont.main <- "CL:0000576"
kass_ref.tbl[kass_ref.tbl$label.fine == "non-classical monocyte",]$label.main <- "monocyte"
kass_ref.tbl[kass_ref.tbl$label.fine == "classical monocyte",]$ont.main <- "CL:0000576"
kass_ref.tbl[kass_ref.tbl$label.fine == "classical monocyte",]$label.main <- "monocyte"
kass_ref.tbl[kass_ref.tbl$label.fine == "plasma cell",]$ont.main <- "CL:0000236"
kass_ref.tbl[kass_ref.tbl$label.fine == "plasma cell",]$label.main <- "B cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "eosinophil",]$ont.main <- "CL:0000771"
kass_ref.tbl[kass_ref.tbl$label.fine == "eosinophil",]$label.main <- "eosinophil"
kass_ref.tbl[kass_ref.tbl$label.fine == "non-classical monocyte",]$ont.main <- "CL:0000576"
kass_ref.tbl[kass_ref.tbl$label.fine == "non-classical monocyte",]$label.main <- "monocyte"
kass_ref.tbl[kass_ref.tbl$label.fine == "naive B cell",]$ont.main <- "CL:0000236"
kass_ref.tbl[kass_ref.tbl$label.fine == "naive B cell",]$label.main <- "B cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "basophil",]$ont.main <- "CL:0000767"
kass_ref.tbl[kass_ref.tbl$label.fine == "basophil",]$label.main <- "basophil"
kass_ref.tbl[kass_ref.tbl$label.fine == "regulatory T cell",]$ont.main <- "CL:0000815"
kass_ref.tbl[kass_ref.tbl$label.fine == "regulatory T cell",]$label.main <- "regulatory T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8-positive, alpha-beta memory T cell",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8-positive, alpha-beta memory T cell",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "fibroblast",]$ont.main <- "CL:0000057"
kass_ref.tbl[kass_ref.tbl$label.fine == "fibroblast",]$label.main <- "fibroblast"
kass_ref.tbl[kass_ref.tbl$label.fine == "class switched memory B cell",]$ont.main <- "CL:0000236"
kass_ref.tbl[kass_ref.tbl$label.fine == "class switched memory B cell",]$label.main <- "B cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "Naive T-helpers",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "Naive T-helpers",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "Central memory T-helpers",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "Central memory T-helpers",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "Transitional memory T-helpers",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "Transitional memory T-helpers",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "T-helpers TEMRA",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "T-helpers TEMRA",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "Effector memory T-helpers",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "Effector memory T-helpers",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD57pos Cytotoxic NK cells",]$ont.main <- "CL:0000623"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD57pos Cytotoxic NK cells",]$label.main <- "natural killer cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD57neg Cytotoxic NK cells",]$ont.main <- "CL:0000623"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD57neg Cytotoxic NK cells",]$label.main <- "natural killer cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD27neg memory B-cells",]$ont.main <- "CL:0000236"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD27neg memory B-cells",]$label.main <- "B cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "Regulatory NK cells",]$ont.main <- "CL:0000623"
kass_ref.tbl[kass_ref.tbl$label.fine == "Regulatory NK cells",]$label.main <- "natural killer cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "Cytotoxic NK cells",]$ont.main <- "CL:0000623"
kass_ref.tbl[kass_ref.tbl$label.fine == "Cytotoxic NK cells",]$label.main <- "natural killer cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "central memory CD8-positive, alpha-beta T cell",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "central memory CD8-positive, alpha-beta T cell",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8+ Ttm",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8+ Ttm",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "effector memory CD8-positive, alpha-beta T cell",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "effector memory CD8-positive, alpha-beta T cell",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "conventional dendritic cell",]$ont.main <- "CL:0000451"
kass_ref.tbl[kass_ref.tbl$label.fine == "conventional dendritic cell",]$label.main <- "dendritic cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "inflammatory macrophage",]$ont.main <- "CL:0000235"
kass_ref.tbl[kass_ref.tbl$label.fine == "inflammatory macrophage",]$label.main <- "macrophage"
kass_ref.tbl[kass_ref.tbl$label.fine == "T follicular helper cell",]$ont.main <- "CL:0000084"
kass_ref.tbl[kass_ref.tbl$label.fine == "T follicular helper cell",]$label.main <- "T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "helper T cell",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "helper T cell",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "endothelial cell",]$ont.main <- "CL:0000115"
kass_ref.tbl[kass_ref.tbl$label.fine == "endothelial cell",]$label.main <- "endothelial cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "T-helper 1 cell",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "T-helper 1 cell",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "T-helper 2 cell",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "T-helper 2 cell",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8+ T-cells PD1 low",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8+ T-cells PD1 low",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8+ T-cells PD1 high",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "CD8+ T-cells PD1 high",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "gamma-delta T cell",]$ont.main <- "CL:0000798"
kass_ref.tbl[kass_ref.tbl$label.fine == "gamma-delta T cell",]$label.main <- "gamma-delta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "alternatively activated macrophage",]$ont.main <- "CL:0000235"
kass_ref.tbl[kass_ref.tbl$label.fine == "alternatively activated macrophage",]$label.main <- "macrophage"
kass_ref.tbl[kass_ref.tbl$label.fine == "macrophage",]$ont.main <- "CL:0000235"
kass_ref.tbl[kass_ref.tbl$label.fine == "macrophage",]$label.main <- "macrophage"
kass_ref.tbl[kass_ref.tbl$label.fine == "T-helper 17 cell",]$ont.main <- "CL:0000912"
kass_ref.tbl[kass_ref.tbl$label.fine == "T-helper 17 cell",]$label.main <- "helper T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "malignant cell",]$ont.main <- "CL:0001064"
kass_ref.tbl[kass_ref.tbl$label.fine == "malignant cell",]$label.main <- "malignant cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "naive thymus-derived CD8-positive, alpha-beta T cell",]$ont.main <- "CL:0000625"
kass_ref.tbl[kass_ref.tbl$label.fine == "naive thymus-derived CD8-positive, alpha-beta T cell",]$label.main <- "CD8-positive, alpha-beta T cell"
kass_ref.tbl[kass_ref.tbl$label.fine == "plasmablast",]$ont.main <- "CL:0000236"
kass_ref.tbl[kass_ref.tbl$label.fine == "plasmablast",]$label.main <- "B cell"


kass_data <- cbind(all_models_expr, laboratory_data_expressions)
kass_data <- kass_data[,kass_ref.tbl$sample]
all(colnames(kass_data) == kass_ref.tbl$sample)


kass.df <- as.data.frame(kass_ref.tbl)
kass.data <- as.matrix(kass_data)
all(kass.df$sample == colnames(kass.data))


saveRDS(kass.df, "/bigdata/almogangel/super_ref_for_xcell2/kass_labels.rds")
saveRDS(kass.data, "/bigdata/almogangel/super_ref_for_xcell2/kass_data.rds")
