library(tidyverse)

celltype_conversion <- read_tsv("celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# petitprez
petitprez.exp <- readRDS("/bigdata/almogangel/xCell2_data/mouse/petitprez/petitprez_tpm.rds")
petitprez.exp <- as.matrix(petitprez.exp)
petitprez.exp <- cbind("Gene" = rownames(petitprez.exp), petitprez.exp)
write.table(petitprez.exp, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/mouse/petitprez_expressions.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

petitprez.truth <- readRDS("/bigdata/almogangel/xCell2_data/mouse/petitprez/petitprez_facs.rds")
celltype <- rownames(petitprez.truth)
celltype <- plyr::mapvalues(celltype, celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
celltype[celltype  == "MAST cells"] <- "mast cell"
petitprez.truth <- as.matrix(apply(petitprez.truth, 2, as.numeric))
rownames(petitprez.truth) <- celltype
petitprez.truth <- petitprez.truth[!rownames(petitprez.truth) %in% c("Monocytes/Macrophages", "NK/T cells", "B derived cells", "B cells germinal"),]
write.table(petitprez.truth, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/mouse/petitprez.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


# load("/bigdata/almogangel/xCell2_data/mouse/dataset_petitprez.rda")
# dataset_petitprez.exp <- as.matrix(dataset_petitprez$expr_mat)
# dataset_petitprez.exp <- cbind("Gene" = rownames(dataset_petitprez.exp), dataset_petitprez.exp)
# write.table(dataset_petitprez.exp, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/mouse/petitprez_expressions.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# dataset_petitprez.truth <- as.matrix(dataset_petitprez$ref)
# dataset_petitprez.truth <- t(dataset_petitprez.truth)
# colnames(dataset_petitprez.truth) <- dataset_petitprez.truth[1,]
# dataset_petitprez.truth <- dataset_petitprez.truth[-1,]
# celltype <- rownames(dataset_petitprez.truth)
# celltype <- plyr::mapvalues(celltype, celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
# dataset_petitprez.truth <- as.matrix(apply(dataset_petitprez.truth, 2, as.numeric))
# rownames(dataset_petitprez.truth) <- celltype
# dataset_petitprez.truth <- dataset_petitprez.truth[rownames(dataset_petitprez.truth) != "Macrophage/Monocyte",]
# write.table(dataset_petitprez.truth, "//bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/mouse/petitprez.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

# chen
chen.exp <- readRDS("/bigdata/almogangel/xCell2_data/mouse/chen/chen_tpm.rds")
chen.exp <- cbind("Gene" = rownames(chen.exp), chen.exp)
write.table(chen.exp, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/mouse/chen_expressions.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


chen.truth <- readRDS("/bigdata/almogangel/xCell2_data/mouse/chen/chen_facs.rds")
rownames(chen.truth) <- plyr::mapvalues(rownames(chen.truth), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
rownames(chen.truth)[1] <- "monocyte"
rownames(chen.truth)[2] <- "B cell"
rownames(chen.truth)[4] <- "CD4-positive, alpha-beta T cell"

t_truth <- colSums(chen.truth[c("CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"),])
chen.truth <- rbind(chen.truth, "T cell" = t_truth)

write.table(chen.truth, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/mouse/chen.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)



