% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xCell2Analysis.R
\name{xCell2Analysis}
\alias{xCell2Analysis}
\title{xCell2Analysis function}
\usage{
xCell2Analysis(
  mix,
  xcell2object,
  minSharedGenes = 0.9,
  rawScores = FALSE,
  spillover = TRUE,
  spilloverAlpha = 0.5,
  BPPARAM = BiocParallel::SerialParam()
)
}
\arguments{
\item{mix}{A bulk mixture of gene expression data (genes in rows, samples in columns).}

\item{xcell2object}{An S4 object of class `xCell2Object`.}

\item{minSharedGenes}{Minimum fraction of shared genes required between the mix and the reference (default: 0.9).}

\item{rawScores}{A Boolean to indicate whether to return raw enrichment scores (default: FALSE).}

\item{spillover}{A Boolean to enable spillover correction on the enrichment scores (default: TRUE). 
Spillover occurs when gene expression patterns overlap between closely related cell types, 
potentially inflating enrichment scores. By correcting for spillover, xCell2 improves the specificity of enrichment scores for closely related cell types. 
`spilloverAlpha` can be tuned to adjust the strength of the correction to avoid under- or over-correction.}

\item{spilloverAlpha}{A numeric value (default: 0.5) that controls the strength of spillover correction. 
Lower values apply weaker correction, while higher values apply stronger correction. An alpha value of 
0.5 is suitable for most cases. However, if your reference contains very closely related cell types, 
consider using a higher value. Conversely, use a lower value if the cell types in your reference are more distinct.}

\item{BPPARAM}{A BiocParallelParam instance that determines the parallelisation strategy. Default is BiocParallel::SerialParam().}
}
\value{
A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
}
\description{
This function performs cell type enrichment analysis to identify proportions of cell types in a bulk gene expression mixture.
}
\examples{
# For detailed example read xCell2 vignette.

# Load "ready to use" xCell2 reference object or generate a new one using `xCell2Train`
data(DICE_demo.xCell2Ref, package = "xCell2")

# Load bulk RNA-Seq gene expression mixture
data(mix_demo, package = "xCell2")

# Run xCell2 cell type enrichment analysis
xcell2_res <- xCell2::xCell2Analysis(mix = mix_demo, xcell2object = DICE_demo.xCell2Ref)

# Example using parallel processing with MulticoreParam
library(BiocParallel)
parallel_param <- MulticoreParam(workers = 2) # Adjust workers as needed
xcell2_res_parallel <- xCell2::xCell2Analysis(
  mix = mix_demo, 
  xcell2object = DICE_demo.xCell2Ref, 
  BPPARAM = parallel_param
)

}
