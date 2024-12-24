#' Subset of the DICE Reference
#'
#' A subset of the DICE reference stored as a `SummarizedExperiment` object for the xCell 2.0 vignette.
#' This demo reference object demonstrates the use of \code{xCell2Train} for generating a custom xCell2 reference.
#'
#' @name dice_demo_ref
#' @docType data
#' @title Subset of the DICE Reference
#' @description Demo reference object derived from the DICE dataset for training xCell2 references.
#' @format A \linkS4class{SummarizedExperiment} object.
#' @usage data(dice_demo_ref, package = "xCell2")
#' @source Schmiedel B et al. (2018).
#' @keywords datasets, references
#' @seealso \code{\link{xCell2Train}} for generating references, and \code{\link{xCell2Analysis}} for enrichment analysis.
NULL

#' Demo Bulk Gene Expression Data (RNA-Seq)
#'
#' A demo mixture matrix for bulk RNA-Seq gene expression data.  
#' Use this dataset to test \code{xCell2Analysis} with pre-trained xCell2 references.
#'
#' @name mix_demo
#' @docType data
#' @title Demo Bulk Gene Expression Data
#' @description Example RNA-Seq data to demonstrate \code{xCell2Analysis}.
#' @format A matrix with genes (rows) and samples (columns).
#' @usage data(mix_demo, package = "xCell2")
#' @keywords datasets, examples
#' @seealso \code{\link{xCell2Analysis}} for enrichment analysis.
NULL

#' Demo xCell2 Reference Object from DICE Subset (human)
#'
#' A demo xCell2 reference object derived from a subset of the DICE dataset.
#' Suitable for demonstrating the use of \code{xCell2Analysis}.
#'
#' @name DICE_demo.xCell2Ref
#' @docType data
#' @title Demo xCell2 Reference Object
#' @description Pre-trained xCell2 reference object based on the DICE dataset.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Cell-type-specific gene signatures.}
#'   \item{\code{dependencies}}{Cell type lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(DICE_demo.xCell2Ref, package = "xCell2")
#' @source Schmiedel B et al. (2018).
#' @keywords datasets, references
#' @seealso \code{\link{xCell2Analysis}} for enrichment analysis, and \code{\link{xCell2Train}} for training custom references.
NULL

#' Tumor Microenvironment Compendium Reference (human)
#'
#' An xCell2 reference object created from the Tumor Microenvironment Compendium dataset.
#'
#' @name TMECompendium.xCell2Ref
#' @docType data
#' @title Tumor Microenvironment Compendium Reference
<<<<<<< HEAD
#' @description A reference object for use in `xCell2Analysis`.
#' @format An `xCell2Object` (S4) containing:
=======
#' @description Pre-trained xCell2 reference object for analyzing tumor microenvironments.
#' @format An \code{xCell2Object} containing:
>>>>>>> devel
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Cell-type-specific gene signatures.}
#'   \item{\code{dependencies}}{Cell type lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(TMECompendium.xCell2Ref, package = "xCell2")
#' @details
<<<<<<< HEAD
#' Normalized data used to train this `xCell2` reference object can be found here:
#' https://science.bostongene.com/kassandra/downloads
#' 
#' @source Curated by Zaitsev A (2022) and trained by Angel A, et al. (2024).
#' @references Zaitsev, A., et al. (2022). Precise reconstruction of the TME using bulk RNA-seq and a machine learning algorithm trained on artificial transcriptomes. Cancer Cell, 40(8), 879-894.
#' @keywords datasets
=======
#' Normalized data for training can be accessed at:  
#' \url{https://science.bostongene.com/kassandra/downloads}.
#' @source Curated by Zaitsev A (2022) and trained by Angel A, et al. (2024).
#' @references
#' Zaitsev, A., et al. (2022). Cancer Cell, 40(8), 879-894.
#' @keywords datasets, references
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
>>>>>>> devel
NULL

#' Tabula Muris Blood Reference (mouse)
#'
#' A pre-trained xCell2 reference object based Tabula Muris dataset.
#'
#' @name TabulaMurisBlood.xCell2Ref
#' @docType data
#' @title Tabula Muris Blood Reference
#' @description Pre-trained xCell2 reference for use in \code{xCell2Analysis} or extending via \code{xCell2Train}.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Gene signatures for cell types.}
#'   \item{\code{dependencies}}{Cell lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(TabulaMurisBlood.xCell2Ref, package = "xCell2")
#' @source The Tabula Muris Consortium (2018); trained by Angel A, et al. (2024).
#' @keywords datasets
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
NULL

#' Tabula Sapiens Blood Reference (human)
#'
#' A pre-trained xCell2 reference object based on the Tabula Sapiens dataset.
#'
#' @name TabulaSapiensBlood.xCell2Ref
#' @docType data
#' @title Tabula Sapiens Blood Reference
#' @description Pre-trained xCell2 reference for use in \code{xCell2Analysis} or extending via \code{xCell2Train}.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Gene signatures for cell types.}
#'   \item{\code{dependencies}}{Cell lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(TabulaSapiensBlood.xCell2Ref, package = "xCell2")
#' @source The Tabula Sapiens Consortium (2022); trained by Angel A, et al. (2024).
#' @keywords datasets
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
NULL

#' Mouse RNA-Seq Data Reference
#'
#' A pre-trained xCell2 reference object based on the MouseRNAseqData dataset.
#'
#' @name MouseRNAseqData.xCell2Ref
#' @docType data
#' @title Mouse RNA-Seq Data Reference
#' @description Pre-trained xCell2 reference for use in \code{xCell2Analysis} or extending via \code{xCell2Train}.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Gene signatures for cell types.}
#'   \item{\code{dependencies}}{Cell lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(MouseRNAseqData.xCell2Ref, package = "xCell2")
#' @source Benayoun B (2019); trained by Angel A, et al. (2024).
#' @keywords datasets
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
NULL

#' PanCancer Reference (human)
#'
#' A pre-trained xCell2 reference object based on the PanCancer dataset for cancer-specific analyses.
#'
#' @name PanCancer.xCell2Ref
#' @docType data
#' @title PanCancer Reference
#' @description Pre-trained xCell2 reference for use in \code{xCell2Analysis} or extending via \code{xCell2Train}.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Gene signatures for cell types.}
#'   \item{\code{dependencies}}{Cell lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(PanCancer.xCell2Ref, package = "xCell2")
#' @source Nofech-Mozes I (2023); trained by Angel A, et al. (2024).
#' @keywords datasets
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
NULL

#' LM22 Reference (human)
#'
#' A pre-trained xCell2 reference object based on the LM22 dataset.
#'
#' @name LM22.xCell2Ref
#' @docType data
#' @title LM22 Reference
#' @description Pre-trained xCell2 reference for use in \code{xCell2Analysis} or extending via \code{xCell2Train}.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Gene signatures for cell types.}
#'   \item{\code{dependencies}}{Cell lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(LM22.xCell2Ref, package = "xCell2")
#' @source Newman AM (2015); trained by Angel A, et al. (2024).
#' @keywords datasets
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
NULL

#' Immune Compendium Reference (human)
#'
#' A pre-trained xCell2 reference object based on the Immune Compendium dataset for immune cell profiling.
#'
#' @name ImmuneCompendium.xCell2Ref
#' @docType data
#' @title Immune Compendium Reference
#' @description Pre-trained xCell2 reference for use in \code{xCell2Analysis} or extending via \code{xCell2Train}.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Gene signatures for cell types.}
#'   \item{\code{dependencies}}{Cell lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(ImmuneCompendium.xCell2Ref, package = "xCell2")
#' @source Curated by Zaitsev A (2022); trained by Angel A, et al. (2024).
#' @keywords datasets
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
NULL

#' Immunologic Genome Project Reference
#'
#' A pre-trained xCell2 reference object based on the Immunologic Genome Project dataset.
#'
#' @name ImmGenData.xCell2Ref
#' @docType data
#' @title Immunologic Genome Project Reference
#' @description Pre-trained xCell2 reference for use in \code{xCell2Analysis} or extending via \code{xCell2Train}.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Gene signatures for cell types.}
#'   \item{\code{dependencies}}{Cell lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(ImmGenData.xCell2Ref, package = "xCell2")
#' @source The Immunological Genome Project Consortium (2008), curated by Aran D (2019); trained by Angel A, et al. (2024).
#' @keywords datasets
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
NULL

#' Blueprint and ENCODE Projects Reference (human)
#'
#' A pre-trained xCell2 reference object based on the Blueprint and ENCODE projects datasets.
#'
#' @name BlueprintEncode.xCell2Ref
#' @docType data
#' @title Blueprint and ENCODE Projects Reference
#' @description Pre-trained xCell2 reference for use in \code{xCell2Analysis} or extending via \code{xCell2Train}.
#' @format An \code{xCell2Object} containing:
#' \describe{
#'   \item{\code{params}}{Linear transformation parameters.}
#'   \item{\code{signatures}}{Gene signatures for cell types.}
#'   \item{\code{dependencies}}{Cell lineage dependencies.}
#'   \item{\code{spill_mat}}{Spillover correction matrix.}
#'   \item{\code{genes_used}}{Genes included in the reference.}
#' }
#' @usage data(BlueprintEncode.xCell2Ref, package = "xCell2")
#' @source Martens JHA and Stunnenberg HG (2013); The ENCODE Project Consortium (2012), curated by Aran D (2019); trained by Angel A, et al. (2024).
#' @keywords datasets
#' @seealso \code{\link{xCell2Analysis}} and \code{\link{xCell2Train}}.
NULL
