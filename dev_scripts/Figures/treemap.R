library(tidyverse)
library(RColorBrewer)
library(treemapify)
library(gridExtra)


refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")


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
                "myeloid suppressor cell",
                "plasmacytoid dendritic cell, human",
                "professional antigen presenting cell",
                "conventional dendritic cell",
                "dendritic cell")

stroma_cells <- c("endothelial cell",
                  "fibroblast",
                  "epithelial cell",
                  "malignant cell")


nk_cells <- c("natural killer cell",
              "mature NK T cell",
              "CD57neg Cytotoxic NK cells",
              "CD57pos Cytotoxic NK cells",
              "conventional dendritic cell",
              "Cytotoxic NK cells",
              "Regulatory NK cells")




celltypes <- sort(table(unlist(refval.tbl$shared_celltypes)), decreasing = TRUE)
t_cells <- t_cells[t_cells %in% names(celltypes)]
b_cells <- b_cells[b_cells %in% names(celltypes)]
myelocytes <- myelocytes[myelocytes %in% names(celltypes)]
stroma_cells <- stroma_cells[stroma_cells %in% names(celltypes)]
nk_cells <- nk_cells[nk_cells %in% names(celltypes)]

all(c(t_cells, b_cells, myelocytes, stroma_cells, nk_cells) %in% names(celltypes))





t.data <- as.data.frame(celltypes[t_cells])
colnames(t.data) <- c("cell_type", "count")
t.data$cell_type <- as.character(t.data$cell_type)
others <- which(t.data$count == 1)
t.data <- t.data[-others,]
t.data <- rbind(t.data, data.frame("cell_type" = "Other T sub-types", "count" = length(others)))
t.data$cell_type <- paste0(t.data$cell_type, " (", t.data$count, ")")
t.data$cell_type <- as.factor(t.data$cell_type)
rownames(t.data) <- t.data$cell_type
t.data$type <- "T cells"

b.data <- as.data.frame(celltypes[b_cells])
colnames(b.data) <- c("cell_type", "count")
b.data$cell_type <- as.character(b.data$cell_type)
b.data$cell_type <- paste0(b.data$cell_type, " (", b.data$count, ")")
b.data$cell_type <- as.factor(b.data$cell_type)
rownames(b.data) <- b.data$cell_type
b.data$type <- "B cells"


myelo.data <- as.data.frame(celltypes[myelocytes])
colnames(myelo.data) <- c("cell_type", "count")
myelo.data$cell_type <- as.character(myelo.data$cell_type)
myelo.data$cell_type <- paste0(myelo.data$cell_type, " (", myelo.data$count, ")")
myelo.data$cell_type <- as.factor(myelo.data$cell_type)
rownames(myelo.data) <- myelo.data$cell_type
myelo.data$type <- "Myelocytes"


stroma.data <- as.data.frame(celltypes[stroma_cells])
colnames(stroma.data) <- c("cell_type", "count")
stroma.data$cell_type <- as.character(stroma.data$cell_type)
stroma.data$cell_type <- paste0(stroma.data$cell_type, " (", stroma.data$count, ")")
stroma.data$cell_type <- as.factor(stroma.data$cell_type)
rownames(stroma.data) <- stroma.data$cell_type
stroma.data$type <- "Stroma"

nk.data <- as.data.frame(celltypes[nk_cells])
colnames(nk.data) <- c("cell_type", "count")
nk.data$cell_type <- as.character(nk.data$cell_type)
others <- which(nk.data$cell_type %in% c("Cytotoxic NK cells", "CD57neg Cytotoxic NK cells", "CD57pos Cytotoxic NK cells"))
nk.data <- nk.data[-others,]
nk.data <- rbind(nk.data, data.frame("cell_type" = "Cytotoxic NK cells", "count" = length(others)))
nk.data$cell_type <- paste0(nk.data$cell_type, " (", nk.data$count, ")")
nk.data$cell_type <- as.factor(nk.data$cell_type)
rownames(nk.data) <- nk.data$cell_type
nk.data$type <- "NK cells"



data.combined <- rbind(t.data, b.data, myelo.data, stroma.data, nk.data)



# Define the base colors for each "type"
type_base_colors <- brewer.pal(n = length(unique(data.combined$type)), name = "Set3")
names(type_base_colors) <- unique(data.combined$type)

data.combined$fill_color <- NA
for (x in (unique(data.combined$type))) {
  type.index <- data.combined$type == x
  data.combined[type.index,]$fill_color <- colorRampPalette(c(type_base_colors[x], "black"))(sum(type.index)+10)[1:sum(type.index)]
}

name_mapping <- c(
  "CD8-positive, alpha-beta T cell (68)" = "CD8+\n(68)",
  "CD4-positive, alpha-beta T cell (52)" = "CD4+\n(52)",
  "T cell (24)" = "T cell\ngeneral (24)",
  "regulatory T cell (32)" = "T-reg (32)",
  "naive thymus-derived CD4-positive, alpha-beta T cell (12)" = "naive CD4+\n(12)",
  "naive thymus-derived CD8-positive, alpha-beta T cell (12)" = "naive CD8+\n(12)",
  "central memory CD8-positive, alpha-beta T cell (10)" = "cm-CD8+\n(10)",
  "effector memory CD8-positive, alpha-beta T cell (10)" = "em-CD8+\n(10)",
  "gamma-delta T cell (5)" = "gd T cell (5)",
  "CD4-positive, alpha-beta memory T cell (4)" = "CD4+ mem\n(4)",
  "CD8-positive, alpha-beta memory T cell (8)" = "CD8+ mem\n(8)",
  "T follicular helper cell (6)" = "Tfh cell (6)",
  "central memory CD4-positive, alpha-beta T cell (2)" = "cm-CD4+\n(2)",
  "effector memory CD4-positive, alpha-beta T cell (2)" = "em-CD4+\n(2)",
  "Other T sub-types (7)" = "Others\n(7)",
  "B cell (28)" = "B cell\n(28)",
  "naive B cell (32)" = "naive B cell\n(32)",
  "memory B cell (21)" = "mem B cell\n(21)",
  "plasma cell (26)" = "plasma cell\n(26)",
  "class switched memory B cell (6)" = "cs-mem B cell\n(6)",
  "Non plasma B-cells (1)" = "non-plasma\nB-cells'\n(1)",
  "plasmablast (3)" = "plasmablast\n(3)",
  "unswitched memory B cell (3)" = "us-mem B cell\n(3)",
  "CD27neg memory B-cells (1)" = "CD27neg\nB cells\n(1)",
  "monocyte (66)" = "monocyte\n(66)",
  "neutrophil (34)" = "neutrophil\n(34)",
  "basophil (4)" = "basophil\n(4)",
  "classical monocyte (10)" = "c-monocyte\n(10)",
  "non-classical monocyte (5)" = "nc-nmonocyte\n(5)",
  "granulocyte (6)" = "granulocyte\n(6)",
  "eosinophil (8)" = "eosinophil\n(8)",
  "macrophage (11)" = "macrophage\n(11)",
  "myeloid cell (3)" = "myeloid cell\n(3)",
  "intermediate monocyte (2)" = "i-monocyte\n(2)",
  "plasmacytoid dendritic cell, human (4)" = "pDC\n(4)",
  "conventional dendritic cell (1)" = "cDC\n(1)",
  "dendritic cell (6)" = "DC\n(6)",
  "endothelial cell (6)" = "endothelial\n(6)",
  "fibroblast (3)" = "fibroblast\n(3)",
  "epithelial cell (2)" = "epithelial\n(2)",
  "malignant cell (4)" = "malignant\n(4)",
  "natural killer cell (40)" = "NK cell\n(40)",
  "conventional dendritic cell (1)1" = "cDC\n(1)",
  "Regulatory NK cells (1)" = "reg-NK\n(1)",
  "Cytotoxic NK cells (3)" = "cyto-NK\n(3)"
)
data.combined$short_cell_type <- name_mapping[rownames(data.combined)]


ggplot(data.combined, aes(area = count, fill = fill_color, label = short_cell_type, subgroup = type)) +
  geom_treemap(color = "white", size = 1) +  # border color and width equivalent
  geom_treemap_text(colour = "black", place = "centre", grow = FALSE, min.size = 6) +  # text customization with bold
  geom_treemap_subgroup_border(colour = "white", size = 2) +  # border for the subgroups
  #geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.5, colour = "black", fontface = "italic", size = 12) +  # text for subgroups
  scale_fill_identity() +  # Use the exact colors defined in fill_color
  guides(fill = "none") +  # remove legend
  ggtitle("") +  # title
  theme(legend.position = "none",  # remove legend
        plot.title = element_text(size = 20, face = "bold"),  # title font size and bold
        aspect.ratio = 0.7)  # aspect ratio equivalent












# Human -------------------
data <- as.data.frame(celltypes[t_cells])
colnames(data) <- c("cell_type", "count")
data$cell_type <- as.character(data$cell_type)
others <- which(data$count == 1)
data <- data[-others,]
data <- rbind(data, data.frame("cell_type" = "Other T sub-types", "count" = length(others)))
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
data$cell_type <- as.factor(data$cell_type)
rownames(data) <- data$cell_type


# Create a mapping from long names to short names
name_mapping <- c(
  "CD8-positive, alpha-beta T cell (68)" = "CD8+\n(68)",
  "CD4-positive, alpha-beta T cell (52)" = "CD4+\n(52)",
  "T cell (24)" = "T cell\ngeneral (24)",
  "regulatory T cell (32)" = "T-reg (32)",
  "naive thymus-derived CD4-positive, alpha-beta T cell (12)" = "naive CD4+\n(12)",
  "naive thymus-derived CD8-positive, alpha-beta T cell (12)" = "naive CD8+\n(12)",
  "central memory CD8-positive, alpha-beta T cell (10)" = "cm-CD8+\n(10)",
  "effector memory CD8-positive, alpha-beta T cell (10)" = "em-CD8+\n(10)",
  "gamma-delta T cell (5)" = "gd T cell (5)",
  "CD4-positive, alpha-beta memory T cell (4)" = "CD4+ mem\n(4)",
  "CD8-positive, alpha-beta memory T cell (8)" = "CD8+ mem\n(8)",
  "T follicular helper cell (6)" = "Tfh cell (6)",
  "central memory CD4-positive, alpha-beta T cell (2)" = "cm-CD4+\n(2)",
  "effector memory CD4-positive, alpha-beta T cell (2)" = "em-CD4+\n(2)",
  "Other T sub-types (7)" = "Others\n(7)"
)
data$short_cell_type <- name_mapping[rownames(data)]

extended_palette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Set3"))(length(unique(data$short_cell_type)))

ggplot(data, aes(area = count, fill = short_cell_type, label = short_cell_type)) +
  geom_treemap(color = "white", size = 1) +  # border color and width equivalent
  geom_treemap_text(colour = "black", place = "centre", reflow = TRUE, size = 8, fontface = "bold") +  # text customization with bold
  scale_fill_manual(values = extended_palette) +  # color palette
  ggtitle("") +  # title
  theme(legend.position = "none",  # remove legend
        plot.title = element_text(size = 20),  # title font size
        aspect.ratio = 1.5)  # aspect ratio equivalent


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
        aspRatio = 1.2)





data <- as.data.frame(celltypes[b_cells])
colnames(data) <- c("cell_type", "count")
data$cell_type <- as.character(data$cell_type)
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
data$cell_type <- as.factor(data$cell_type)

b_tree <- treemap(data,
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
data$cell_type <- as.character(data$cell_type)
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
data$cell_type <- as.factor(data$cell_type)

myelo_tree <- treemap(data,
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
                  aspRatio = 1)


data <- as.data.frame(celltypes[stroma_cells])
colnames(data) <- c("cell_type", "count")
data$cell_type <- as.character(data$cell_type)
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
data$cell_type <- as.factor(data$cell_type)

stroma_tree <- treemap(data,
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


data <- as.data.frame(celltypes[nk_cells])
colnames(data) <- c("cell_type", "count")
data$cell_type <- as.character(data$cell_type)
others <- which(data$cell_type %in% c("Cytotoxic NK cells", "CD57neg Cytotoxic NK cells", "CD57pos Cytotoxic NK cells"))
data <- data[-others,]
data <- rbind(data, data.frame("cell_type" = "Cytotoxic NK cells", "count" = length(others)))
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
data$cell_type <- as.factor(data$cell_type)
others_tree <- treemap(data,
                       index = "cell_type",
                       palette = "Set3",
                       vSize = "count",
                       type = "index",
                       title = "NK cells",
                       fontsize.title = 20,
                       border.col = "white",
                       border.lwds = c(3, 1),
                       fontsize.labels = 10,
                       fontsize.labels.add = 0.5,
                       fontcolor.labels = "black",
                       bg.labels = "transparent",
                       align.labels = list(c("center", "center"), c("center", "center")),
                       aspRatio = 0.5)




# Archive ---------------------

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


all(c(t_cells, b_cells, myelocytes, stroma_cells, nk_cells, others) %in% names(celltypes))




data <- as.data.frame(celltypes[t_cells])
colnames(data) <- c("cell_type", "count")
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
t_tree <- treemap(data,
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
b_tree <- treemap(data,
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
myelo_tree <- treemap(data,
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
stroma_tree <- treemap(data,
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
others_tree <- treemap(data,
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




# Mouse -------------
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val_mouse.rds")
celltypes <- sort(table(unlist(refval.tbl$shared_celltypes)), decreasing = TRUE)


data <- as.data.frame(celltypes)
colnames(data) <- c("cell_type", "count")
data$cell_type <- as.character(data$cell_type)
others <- which(data$count == 1)
data <- data[-others,]
data <- rbind(data, data.frame("cell_type" = "Other T sub-types", "count" = length(others)))
data$cell_type <- paste0(data$cell_type, " (", data$count, ")")
data$cell_type <- as.factor(data$cell_type)

mouse_tree <- treemap(data,
        index = "cell_type",
                  palette = "Set3",
                  vSize = "count",
                  type = "index",
                  title = "",
                  fontsize.title = 20,
                  border.col = "white",
                  border.lwds = c(3, 1),
                  fontsize.labels = 10,
                  fontsize.labels.add = 0.5,
                  fontcolor.labels = "black",
                  bg.labels = "transparent",
                  align.labels = list(c("center", "center"), c("center", "center")),
                  aspRatio = 0.5)



