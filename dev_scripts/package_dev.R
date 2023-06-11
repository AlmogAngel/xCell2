devtools::install_github("AlmogAngel/xCell2", force = TRUE, dep = FALSE)


library(xCell2)

setwd("/bigdata/almogangel/xCell2/")
devtools::document()
devtools::install(pkg = "/bigdata/almogangel/xCell2/")
devtools::check()

?xCell2Analysis()
?xCell2GetLineage()
?xCell2Train()
data("ts_labels_with_ontology")

# data
ts_labels <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/ts_labels.rds")
ts_labels_with_ontology <- ts_labels %>%
  as_tibble()
usethis::use_data(ts_labels_with_ontology)

# license
usethis::use_mit_license()

