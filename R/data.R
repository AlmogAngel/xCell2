#' Human to Mouse Gene Symbol Conversion Table
#'
#' A data frame for converting human gene symbols to mouse gene symbols.
#'
#' @name human_mouse_gene_symbols
#' @docType data
#' @title Human to Mouse Gene Symbols
#' @description A table that provides a mapping of human gene symbols to mouse gene symbols.
#' @format A data frame with 2 columns:
#' \describe{
#'   \item{human}{Character. Human gene symbols.}
#'   \item{mouse}{Character. Corresponding mouse gene symbols.}
#' }
#' @usage data(human_mouse_gene_symbols, package = "xCell2")
#' @keywords datasets
"human_mouse_gene_symbols"

#' A Subset of the DICE Reference
#'
#' A subset `SummarizedExperiment` object of the DICE reference for the xCell 2.0 vignette.
#' This is a demo reference object to learn how to use the `xCell2Train` function.
#'
#' @name dice_demo_ref
#' @docType data
#' @title Subset of the DICE Reference
#' @description A demo reference object derived from the DICE reference for use in the xCell 2.0 vignette.
#' @format A `SummarizedExperiment` object.
#' @usage data(dice_demo_ref, package = "xCell2")
#' @source Schmiedel B et al. (2018)
#' @keywords datasets
NULL

#' Demo Bulk Gene Expression Data (RNA-Seq)
#'
#' A demo mixture matrix for bulk RNA-Seq gene expression data, used to demonstrate how to use the `xCell2Train` function.
#'
#' @name mix_demo
#' @docType data
#' @title Demo Bulk Gene Expression Data
#' @description A demo mixture matrix for bulk RNA-Seq gene expression data. This dataset can be used to learn how to use the `xCell2Train` function.
#' @format A matrix with rows representing genes and columns representing samples.
#' @usage data(mix_demo, package = "xCell2")
#' @keywords datasets
NULL

#' Demo xCell2 Reference Object from a Subset of the DICE Reference
#'
#' A demo xCell2 reference object derived from a subset of the DICE reference dataset.
#'
#' @name DICE_demo.xCell2Ref
#' @docType data
#' @title Demo xCell2 Reference Object
#' @description This object can be used to learn and demonstrate the functionality of the `xCell2Analysis` function.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(DICE_demo.xCell2Ref, package = "xCell2")
#' @source Schmiedel B et al. (2018)
#' @keywords datasets
NULL

#' xCell2 Object Trained from the Tumor Microenvironment Compendium Reference
#'
#' An xCell2 reference object created from the Tumor Microenvironment Compendium reference.
#'
#' @name TMECompendium.xCell2Ref
#' @docType data
#' @title Tumor Microenvironment Compendium Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(TMECompendium.xCell2Ref, package = "xCell2")
#' @source Curated by Zaitsev A (2022) and trained by Angel A, et al. (2024).
#' @keywords datasets
NULL

#' xCell2 Object Trained from the Tabula Muris Reference
#'
#' An xCell2 reference object created from the Tabula Muris reference.
#'
#' @name TabulaMurisBlood.xCell2Ref
#' @docType data
#' @title Tabula Muris Blood Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(TabulaMurisBlood.xCell2Ref, package = "xCell2")
#' @source The Tabula Muris Consortium (2018); trained by Angel A, et al. (2024).
#' @keywords datasets
NULL

#' xCell2 Object Trained from the Tabula Sapiens Reference
#'
#' An xCell2 reference object created from the Tabula Sapiens reference.
#'
#' @name TabulaSapiensBlood.xCell2Ref
#' @docType data
#' @title Tabula Sapiens Blood Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(TabulaSapiensBlood.xCell2Ref, package = "xCell2")
#' @source The Tabula Sapiens Consortium (2022); trained by Angel A, et al. (2024).
#' @keywords datasets
NULL

#' xCell2 Object Trained from the MouseRNAseqData Reference
#'
#' An xCell2 reference object created from the MouseRNAseqData reference.
#'
#' @name MouseRNAseqData.xCell2Ref
#' @docType data
#' @title Mouse RNAseq Data Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(MouseRNAseqData.xCell2Ref, package = "xCell2")
#' @source Generated by Benayoun B (2019); trained by Angel A, et al. (2024).
#' @keywords datasets
NULL

#' xCell2 Object Trained from the PanCancer Reference
#'
#' An xCell2 reference object created from the PanCancer reference.
#'
#' @name PanCancer.xCell2Ref
#' @docType data
#' @title PanCancer Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(PanCancer.xCell2Ref, package = "xCell2")
#' @source Generated by Nofech-Mozes I (2023); trained by Angel A, et al. (2024).
#' @keywords datasets
NULL

#' xCell2 Object Trained from the LM22 Reference
#'
#' An xCell2 reference object created from the LM22 reference.
#'
#' @name LM22.xCell2Ref
#' @docType data
#' @title LM22 Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(LM22.xCell2Ref, package = "xCell2")
#' @source Generated by Newman AM (2015); trained by Angel A, et al. (2024).
#' @keywords datasets
NULL

#' xCell2 Object Trained from the Immune Compendium Reference
#'
#' An xCell2 reference object created from the Immune Compendium reference.
#'
#' @name ImmuneCompendium.xCell2Ref
#' @docType data
#' @title Immune Compendium Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(ImmuneCompendium.xCell2Ref, package = "xCell2")
#' @source Curated by Zaitsev A (2022); trained by Angel A, et al. (2024).
#' @keywords datasets
NULL

#' xCell2 Object Trained from the Immunologic Genome Project Reference
#'
#' An xCell2 reference object created from the Immunologic Genome Project reference.
#'
#' @name ImmGenData.xCell2Ref
#' @docType data
#' @title Immunologic Genome Project Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(ImmGenData.xCell2Ref, package = "xCell2")
#' @source The Immunological Genome Project Consortium (2008), curated by Aran D (2019); trained by Angel A, et al. (2024).
#' @keywords datasets
NULL

#' xCell2 Object Trained from the Blueprint and ENCODE Projects Reference
#'
#' An xCell2 reference object created from the Blueprint and ENCODE projects reference.
#'
#' @name BlueprintEncode.xCell2Ref
#' @docType data
#' @title Blueprint and ENCODE Projects Reference
#' @description A reference object for use in `xCell2Analysis` or to create a new reference using `xCell2Train`.
#' @format An `xCell2Object` (S4) containing:
#' \describe{
#'   \item{\code{params}}{A tibble of parameters for linear transformation.}
#'   \item{\code{signatures}}{A list containing cell-type signatures.}
#'   \item{\code{dependencies}}{A list specifying cell type dependencies.}
#'   \item{\code{spill_mat}}{A matrix representing spillover correction factors.}
#'   \item{\code{genes_used}}{A character vector of genes used to generate the reference.}
#' }
#' @usage data(BlueprintEncode.xCell2Ref, package = "xCell2")
#' @source Martens JHA and Stunnenberg HG (2013); The ENCODE Project Consortium (2012), curated by Aran D (2019); trained by Angel A, et al. (2024).
#' @keywords datasets
NULL
