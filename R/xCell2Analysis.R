#' xCell2Analysis function
#'
#' This function performs cell type enrichment analysis to identify proportions of cell types in a bulk gene expression mixture.
#'
#' @importFrom singscore rankGenes simpleScore
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom pracma lsqlincon
#' @param mix A bulk mixture of gene expression data (genes in rows, samples in columns).
#' @param xcell2object An S4 object of class `xCell2Object`.
#' @param minSharedGenes Minimum fraction of shared genes required between the mix and the reference (default: 0.9).
#' @param rawScores Boolean to indicate whether to return raw enrichment scores (default: FALSE).
#' @param spillover Boolean - use spillover correction on the transformed enrichment scores? (default: TRUE).
#' @param spilloverAlpha A numeric for spillover alpha value (default: 0.5).
#' @param numThreads Number of threads for parallel processing (default: 1).
#' @examples
#' # For detailed example read xCell2 vignette.
#' 
#' # Load "ready to use" xCell2 reference object or generate a new one using `xCell2Train`
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' 
#' # Load bulk RNA-Seq gene expression mixture
#' data(mix_demo, package = "xCell2")
#' 
#' # Run xCell2 cell type enrichment analysis
#' xcell2_res <- xCell2::xCell2Analysis(mix = mix_demo, xcell2object = DICE_demo.xCell2Ref)
#' 
#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2.
#' @export
xCell2Analysis <- function(mix,
                           xcell2object,
                           minSharedGenes = 0.9,
                           rawScores = FALSE,
                           spillover = TRUE,
                           spilloverAlpha = 0.5,
                           numThreads = 1) {

  scoreMix <- function(cellType, mixRanked, signaturesCellType) {
    scores <- vapply(signaturesCellType, function(sig) {
      singscore::simpleScore(mixRanked, upSet = sig, centerScore = FALSE)$TotalScore
    }, FUN.VALUE = double(ncol(mixRanked)))
    rownames(scores) <- colnames(mixRanked)
    return(scores)
  }

  param <- BiocParallel::MulticoreParam(workers = numThreads)

  # Check reference/mixture genes intersection
  genesIntersectFrac <- length(intersect(rownames(mix), getGenesUsed(xcell2object))) /
    length(getGenesUsed(xcell2object))
  if (genesIntersectFrac < minSharedGenes) {
    stop("This xCell2 reference shares ",
         genesIntersectFrac, " genes with the mixtures and minSharedGenes = ",
         minSharedGenes, ".",
         "\n", "Consider training a new xCell2 reference or adjusting minSharedGenes.")
  }

  # Rank mix gene expression matrix
  mixRanked <- singscore::rankGenes(mix[getGenesUsed(xcell2object), ])

  # Score and predict
  sigsCellTypes <- unique(unlist(lapply(names(getSignatures(xcell2object)),
                                        function(x) { strsplit(x, "#")[[1]][1] })))

  # Get raw enrichment scores
  resRaw <- BiocParallel::bplapply(sigsCellTypes, function(cellType) {
    signaturesCellType <- getSignatures(xcell2object)[startsWith(names(
      getSignatures(xcell2object)), paste0(cellType, "#"))]
    scores <- scoreMix(cellType, mixRanked, signaturesCellType)
    return(scores)
  }, BPPARAM = param)

  names(resRaw) <- sigsCellTypes

  res <- t(vapply(resRaw, function(cellTypeScores) {
    rowMeans(cellTypeScores)
  }, FUN.VALUE = double(nrow(resRaw[[1]]))))

  if (rawScores) {
    return(res)
  } else {
    # Linear transformation
    res <- t(vapply(rownames(res), function(cellType) {
      cellTypeRes <- res[cellType, ]

      # Linear transformation
      a <- getParams(xcell2object)[xcell2object@params$celltype == cellType,]$a
      b <- getParams(xcell2object)[xcell2object@params$celltype == cellType,]$b
      m <- getParams(xcell2object)[xcell2object@params$celltype == cellType,]$m

      cellTypeRes <- (cellTypeRes^(1/b)) / a
      cellTypeRes <- cellTypeRes * m

      # Shift values
      cellTypeRes <- cellTypeRes - min(cellTypeRes)
      cellTypeRes <- round(cellTypeRes, 5)

      return(cellTypeRes)
    }, FUN.VALUE = double(ncol(res))))
  }

  if (spillover) {
    # Spillover correction
    spillMat <- getSpillMat(xcell2object) * spilloverAlpha
    diag(spillMat) <- 1

    rows <- intersect(rownames(res), rownames(spillMat))

    scoresCorrected <- apply(res[rows, ], 2, function(x) pracma::lsqlincon(spillMat[rows, rows],
                                                                           x, lb = 0))
    scoresCorrected[scoresCorrected < 0] <- 0
    rownames(scoresCorrected) <- rows

    return(scoresCorrected)
  } else {
    return(res)
  }
}
