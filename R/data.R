#' @name hs.genelist
#' @title List of human genes that can be removed from ref
#' @description Genes include 7 groups: "Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX" and "chrY" (Credit from GitHub: BayesPrism/inst/extdata/genelist.hs.new.txt)
#' @format A data frame with 3 columns:
#' \describe{
#'   \item{gene_group}{gene groups}
#'   \item{ensembl}{gene Ensemble ID}
#'   \item{symbol}{gene symbol}
#' }
"hs.genelist"

#' @name celltype.data
#' @title aa
#' @description Gaaa
#' @format A tibble with 4 columns:
#' \describe{
#'   \item{xCell2_labels}{xCell2_labels}
#'   \item{all_labels}{all_labels}
#'   \item{ont}{ontl}
#'   \item{essential_genes}{essential_genes}
#' }
"celltype.data"
