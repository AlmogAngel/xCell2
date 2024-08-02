# Load necessary library
library(treemap)

cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

celltypes <- sort(table(c(unlist(lapply(cyto.vals$truth$tumor, rownames)),
               unlist(lapply(cyto.vals$truth$other, rownames)),
               unlist(lapply(cyto.vals$truth$blood, rownames)))), decreasing = T)

celltypes <- celltypes[!names(celltypes) %in% c("Unclassified", "Unidentified", "Other cells",
                                                "B Ex+SM", "B Naive+Memory", "B Naive+NSM",
                                                "Monocytes C+I", "Monocytes NC+I", "T CD4 Naive+Tregs",
                                                "T CD8 EM+TE", "Tfh+Th", "Tfh+Th1-17", "CD44+Unidentified",
                                                "Monocytes.NC+I", "Other_T_helpers", "Other", "Innate", "T Innate",
                                                "Adaptive", "General_cells", "Myeloid Phagocytes", "Progenitors", "T Innate",
                                                "Th1/Th17", "Th17_Th22", "Neutrophils LD", "Non plasmablasts", "Th1_Th17")]


t_cells <- c("CD8-positive, alpha-beta T cell",
             "CD4-positive, alpha-beta T cell",
             "T cell",
             "regulatory T cell",
             "naive thymus-derived CD4-positive, alpha-beta T cell",
             "naive thymus-derived CD8-positive, alpha-beta T cell",
             "central memory CD8-positive, alpha-beta T cell",
             "effector memory CD8-positive, alpha-beta T cell",
             "gamma-delta T cell",
             "CD4-positive, alpha-beta memory T cell",
             "CD8-positive, alpha-beta memory T cell",
             "helper T cell",
             "memory T cell",
             "naive T cell",
             "T follicular helper cell",
             "T-helper 1 cell",
             "T-helper 2 cell",
             "CD4+ Tte",
             "CD8+ T-cells PD1 high",
             "CD8+ T-cells PD1 low",
             "CD8+ Tte",
             "central memory CD4-positive, alpha-beta T cell",
             "effector CD4-positive, alpha-beta T cell",
             "effector CD8-positive, alpha-beta T cell",
             "effector memory CD4-positive, alpha-beta T cell",
             "mucosal invariant T cell",
             "T gd non-Vd2",
             "T gd Vd2",
             "T-helper 17 cell",
             "CD8+ Ttm",
             "Central memory T-helpers",
             "Effector memory T-helpers",
             "Memory T-helpers",
             "Naive T-helpers",
             "T-helpers TEMRA",
             "Transitional memory T-helpers")

b_cells <- c("B cell",
             "naive B cell",
             "memory B cell",
             "plasma cell",
             "B cells",
             "class switched memory B cell",
             "Non plasma B-cells",
             "plasmablast",
             "unswitched memory B cell",
             "pro-B cell",
             "CD27neg memory B-cells",
             "exhausted B cell")


myelocytes <- c("monocyte",
                "neutrophil",
                "basophil",
                "classical monocyte",
                "non-classical monocyte",
                "granulocyte",
                "eosinophil",
                "macrophage",
                "myeloid cell",
                "myeloid dendritic cell",
                "intermediate monocyte",
                "phagocyte",
                "myeloid suppressor cell")

stroma_cells <- c("endothelial cell",
                  "fibroblast",
                  "epithelial cell")


# nk_cells <- c("natural killer cell",
#               "mature NK T cell",
#               "CD57neg Cytotoxic NK cells",
#               "CD57pos Cytotoxic NK cells",
#               "conventional dendritic cell",
#               "Cytotoxic NK cells",
#               "Regulatory NK cells")

others <- c("lymphocyte",
            "dendritic cell",
            "malignant cell",
            "plasmacytoid dendritic cell, human",
            "professional antigen presenting cell",
            "natural killer cell",
            "mature NK T cell",
            "CD57neg Cytotoxic NK cells",
            "CD57pos Cytotoxic NK cells",
            "conventional dendritic cell",
            "Cytotoxic NK cells",
            "Regulatory NK cells")


all(c(t_cells, b_cells, myelocytes, stroma_cells, nk_cells, others) %in% names(celltypes))




data <- as.data.frame(celltypes[t_cells])
colnames(data) <- c("cell_type", "count")
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
treemap(data,
        index = "cell_type",
        palette = "Set3",
        vSize = "count",
        type = "index",
        title = "T cells",
        fontsize.title = 20,
        border.col = "white",
        border.lwds = c(3, 1),
        fontsize.labels = 10,
        fontsize.labels.add = 0.5,
        fontcolor.labels = "black",
        bg.labels = "transparent",
        align.labels = list(c("center", "center"), c("center", "center")),
        aspRatio = 1.3)

data <- as.data.frame(celltypes[b_cells])
colnames(data) <- c("cell_type", "count")
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
treemap(data,
        index = "cell_type",
        palette = "Set3",
        vSize = "count",
        type = "index",
        title = "B cells",
        fontsize.title = 20,
        border.col = "white",
        border.lwds = c(3, 1),
        fontsize.labels = 10,
        fontsize.labels.add = 0.5,
        fontcolor.labels = "black",
        bg.labels = "transparent",
        align.labels = list(c("center", "center"), c("center", "center")),
        aspRatio = 0.9)

data <- as.data.frame(celltypes[myelocytes])
colnames(data) <- c("cell_type", "count")
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
treemap(data,
        index = "cell_type",
        palette = "Set3",
        vSize = "count",
        type = "index",
        title = "Myeloids",
        fontsize.title = 20,
        border.col = "white",
        border.lwds = c(3, 1),
        fontsize.labels = 10,
        fontsize.labels.add = 0.5,
        fontcolor.labels = "black",
        bg.labels = "transparent",
        align.labels = list(c("center", "center"), c("center", "center")),
        aspRatio = 0.9)

data <- as.data.frame(celltypes[stroma_cells])
colnames(data) <- c("cell_type", "count")
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
treemap(data,
        index = "cell_type",
        palette = "Set3",
        vSize = "count",
        type = "index",
        title = "Stroma",
        fontsize.title = 20,
        border.col = "white",
        border.lwds = c(3, 1),
        fontsize.labels = 10,
        fontsize.labels.add = 0.5,
        fontcolor.labels = "black",
        bg.labels = "transparent",
        align.labels = list(c("center", "center"), c("center", "center")),
        aspRatio = 0.5)


data <- as.data.frame(celltypes[others])
colnames(data) <- c("cell_type", "count")
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
treemap(data,
        index = "cell_type",
        palette = "Set3",
        vSize = "count",
        type = "index",
        title = "Other",
        fontsize.title = 20,
        border.col = "white",
        border.lwds = c(3, 1),
        fontsize.labels = 10,
        fontsize.labels.add = 0.5,
        fontcolor.labels = "black",
        bg.labels = "transparent",
        align.labels = list(c("center", "center"), c("center", "center")),
        aspRatio = 0.75)












