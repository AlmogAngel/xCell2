devtools::install_github("AlmogAngel/xCell2")
library(xCell2)

xCell2Analysis()
data("ts_labels_with_ontology")


# data
ts_labels <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/ts_labels.rds")
ts_labels_with_ontology <- ts_labels %>%
  as_tibble()
usethis::use_data(ts_labels_with_ontology)

# dependencies
usethis::use_package("Seurat")
usethis::use_package("kneedle")
usethis::use_package("Rfast")
usethis::use_package("sparseMatrixStats")
usethis::use_package("GSEABase")
usethis::use_package("singscore")
usethis::use_package("ontoProc")
usethis::use_package("ontologyIndex")
usethis::use_package("tibble")
usethis::use_package("dplyr")
usethis::use_package("plyr")

devtools::document()
