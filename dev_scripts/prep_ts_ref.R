library(tidyverse)
source("R/xCell2.R")

ts <- readRDS("/bigdata/almogangel/TabulaSapiens/tabula_sapiens.rds")
ts@assays$RNA@key <- "rna_"

# Use only 10X data
ts <- subset(ts, subset = assay == "10x 3' v3")

# Make labels
ts_labels <- tibble(ont = ts@meta.data$cell_type_ontology_term_id, label = ts@meta.data$cell_type,
                    sample = rownames(ts@meta.data), dataset = ts@meta.data$donor, tissue = ts@meta.data$tissue)

# Make main labels and ontology using sref
sref <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/human_bulk_ref.rds")
sref_labels.tbl <- as_tibble(sref@labels)
rm(sref)

ts_labels <- ts_labels %>%
  rowwise() %>%
  rename(label.fine = label,
         ont.fine = ont) %>%
  mutate(label.main = plyr::mapvalues(label.fine, from = sref_labels.tbl$label.fine, to = sref_labels.tbl$label.main, warn_missing = FALSE),
         ont.main = plyr::mapvalues(ont.fine, from = sref_labels.tbl$ont.fine, to = sref_labels.tbl$ont.main, warn_missing = FALSE)) %>%
  select(ont.main, label.main, ont.fine, label.fine, sample, dataset, tissue) %>%
  mutate_all(as.character)

ts_labels.new <- ts_labels

# Manually fix main labels and ontology
ct2change <- as.character(unique(ts_labels.new$label.main[!ts_labels.new$label.main %in% sref_labels.tbl$label.main]))

ts_labels.new[ts_labels.new$label.fine == "T follicular helper cell",]$ont.main <- "CL:0002038"
ts_labels.new[ts_labels.new$label.fine == "T follicular helper cell",]$label.main <- "T follicular helper cell"

x <- ct2change[grep("dendritic cell", ct2change)]
x <- c(x, "Langerhans cell")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000451"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "dendritic cell"
x <- ct2change[grep("NK T cell", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000814"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "mature NK T cell"
x <- ct2change[grep("cholangiocyte", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:1000488"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "cholangiocyte"
x <- c("mucus secreting cell", "serous cell of epithelium of trachea", "serous cell of epithelium of bronchus")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000151"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "secretory cell"
x <- ct2change[grep("CD4", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000624"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "CD4-positive, alpha-beta T cell"
x <- ct2change[grep("CD8", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000625"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "CD8-positive, alpha-beta T cell"
x <- ct2change[grep("goblet", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000160"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "goblet cell"
x <- ct2change[grep("endothelial", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000115"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "endothelial cell"
x <- ct2change[grep("epithelial", ct2change)]
x <- x[x != "ciliated epithelial cell" ]
x <- c(x, "luminal cell of prostate epithelium")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0002632"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "epithelial cell"
x <- ct2change[grep("ciliated", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000064"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "ciliated cell"
x <- ct2change[grep("muscle", ct2change)]
x <- x[x != "skeletal muscle satellite stem cell" ]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000187"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "muscle cell"
x <- c("fibroblast of breast")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000057"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "fibroblast"
x <- c("retinal bipolar neuron")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000540"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "neuron"
x <- c("radial glial cell", "Muller cell")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000125"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "glial cell"
x <- ct2change[grep("stem cell", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000034"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "stem cell"
x <- ct2change[grep("basal cell", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000646"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "basal cell"
x <- ct2change[grep("acinar cell", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000622"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "acinar cell"
x <- c("naive regulatory T cell")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000815"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "regulatory T cell"
x <- c("enterocyte of epithelium of large intestine", "enterocyte of epithelium of small intestine")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000584"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "enterocyte"
x <- c("paneth cell of colon", "paneth cell of epithelium of small intestine")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000510"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "paneth cell"
x <- c("intestinal enteroendocrine cell")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000164"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "enteroendocrine cell"
x <- c("immature natural killer cell")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000623"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "natural killer cell"
x <- ct2change[grep("ionocyte", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0005006"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "ionocyte"
x <- ct2change[grep("thymocyte", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000893"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "thymocyte"
x <- c("intestinal tuft cell")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0002204"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "brush cell"
x <- ct2change[grep("pneumocyte", ct2change)]
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000322"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "pneumocyte"
x <- c("retina horizontal cell", "retinal ganglion cell")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000540"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "neuron"
x <- c("duodenum glandular cell", "enteroendocrine", "acinar cell")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0000150"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "glandular epithelial cell"
x <- c("transit amplifying cell of colon", "transit amplifying cell of small intestine")
ts_labels.new[ts_labels.new$label.main %in% x,]$ont.main <- "CL:0009010"
ts_labels.new[ts_labels.new$label.main %in% x,]$label.main <- "transit amplifying cell"

# All tissues -----
# Make ref
ts_ref <- ts@assays$RNA@counts[,ts_labels.new$sample]
# Use symbols as gene IDs
rownames(ts_ref) <- ts@assays$RNA@meta.features$feature_name

ts_labels <- as.data.frame(ts_labels.new)
all(ts_labels$sample == colnames(ts_ref))

saveRDS(ts_labels, "/bigdata/almogangel/super_ref_for_xcell2/ts_labels.rds")
saveRDS(ts_ref, "/bigdata/almogangel/super_ref_for_xcell2/ts_data.rds")

# Blood -----

ts_labels_blood <- ts_labels.new %>%
  filter(tissue == "blood")

# Manually add cell types that are in the blood validation datasets

table(ts_labels_blood[ts_labels_blood$label.main == "B cell",]$label.fine)
ts_labels_blood[ts_labels_blood$label.main == "B cell",]$label.main <- ts_labels_blood[ts_labels_blood$label.main == "B cell",]$label.fine
ts_labels_blood[ts_labels_blood$label.main == "B cell",]$ont.main <- ts_labels_blood[ts_labels_blood$label.main == "B cell",]$ont.fine


table(ts_labels_blood[ts_labels_blood$label.main == "mature NK T cell",]$label.fine)
ts_labels_blood[ts_labels_blood$label.main == "mature NK T cell",]$label.main <- ts_labels_blood[ts_labels_blood$label.main == "mature NK T cell",]$label.fine
ts_labels_blood[ts_labels_blood$label.main == "mature NK T cell",]$ont.main <- ts_labels_blood[ts_labels_blood$label.main == "mature NK T cell",]$ont.fine

table(ts_labels_blood[ts_labels_blood$label.main == "monocyte",]$label.fine)
ts_labels_blood[ts_labels_blood$label.main == "monocyte",]$label.main <- ts_labels_blood[ts_labels_blood$label.main == "monocyte",]$label.fine
ts_labels_blood[ts_labels_blood$label.main == "monocyte",]$ont.main <- ts_labels_blood[ts_labels_blood$label.main == "monocyte",]$ont.fine

table(ts_labels_blood[ts_labels_blood$label.main == "CD4-positive, alpha-beta T cell",]$label.fine)
ts_labels_blood[ts_labels_blood$label.main == "CD4-positive, alpha-beta T cell",]$label.main <- ts_labels_blood[ts_labels_blood$label.main == "CD4-positive, alpha-beta T cell",]$label.fine
ts_labels_blood[ts_labels_blood$label.main == "CD4-positive, alpha-beta T cell",]$ont.main <- ts_labels_blood[ts_labels_blood$label.main == "CD4-positive, alpha-beta T cell",]$ont.fine

table(ts_labels_blood[ts_labels_blood$label.main == "CD8-positive, alpha-beta T cell",]$label.fine)
ts_labels_blood[ts_labels_blood$label.main == "CD8-positive, alpha-beta T cell",]$label.main <- ts_labels_blood[ts_labels_blood$label.main == "CD8-positive, alpha-beta T cell",]$label.fine
ts_labels_blood[ts_labels_blood$label.main == "CD8-positive, alpha-beta T cell",]$ont.main <- ts_labels_blood[ts_labels_blood$label.main == "CD8-positive, alpha-beta T cell",]$ont.fine


# Make ref
ts_ref <- ts@assays$RNA@counts[,ts_labels_blood$sample]
# Use symbols as gene IDs
rownames(ts_ref) <- ts@assays$RNA@meta.features$feature_name

ts_labels <- as.data.frame(ts_labels_blood)
all(ts_labels$sample == colnames(ts_ref))

saveRDS(ts_labels, "/bigdata/almogangel/super_ref_for_xcell2/ts_bloodlabels.rds")
saveRDS(ts_ref, "/bigdata/almogangel/super_ref_for_xcell2/ts_blood_data.rds")
