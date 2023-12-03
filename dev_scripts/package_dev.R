devtools::install_github("AlmogAngel/xCell2", force = TRUE, dep = FALSE)

library(xCell2)
library(tidyverse)

setwd("/bigdata/almogangel/xCell2/")
devtools::document()
devtools::check()


?xCell2Analysis()
?xCell2GetLineage()
?xCell2Train()
data("ts_labels_with_ontology")

# data
ts_labels <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/ts_labels.rds")
ts_labels_with_ontology <- ts_labels %>%
  as_tibble()

gencode.v22.broad.category <- read_tsv("/bigdata/almogangel/xCell2/data/gencode.v22.broad.category.txt", col_names = F)
colnames(gencode.v22.broad.category) <- c("chrm", "start", "end", "ensembl1", "symbol", "strand", "info", "ensembl", "genes_type")
gencode.v22.broad.category <- gencode.v22.broad.category %>%
  select(genes_type, ensembl, symbol) %>%
  as.data.frame()

celltype.data <- read_tsv("/bigdata/almogangel/xCell2_data/celltypes.data.tsv") %>%
  rowwise() %>%
  mutate(all_labels = paste0(xCell2_labels, ";", all_labels)) %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  mutate(essential_genes = str_split(essential_genes, ";")) %>%
  unnest(cols = c(all_labels))

usethis::use_data(ts_labels_with_ontology)
usethis::use_data(hs.genelist)
usethis::use_data(celltype.data, overwrite = T)

usethis::use_data(gencode.v22.broad.category, overwrite = T)


# license
usethis::use_mit_license()
