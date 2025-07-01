#' Perform Cell Type Enrichment Analysis
#'
#' Estimates the relative enrichment of cell types in a bulk gene expression mixture.
#' This function uses gene signatures from a pre-trained \code{xCell2Object} to compute enrichment scores, 
#' with options for linear transformation and spillover correction to improve specificity.
#'
#' @importFrom magrittr %>%
#' @importFrom singscore rankGenes simpleScore getStableGenes
#' @importFrom BiocParallel MulticoreParam SerialParam bplapply
#' @importFrom pracma lsqlincon
#' @importFrom progress progress_bar
#' @importFrom quadprog solve.QP
#' 
#' @param mix A bulk mixture of gene expression matrix (genes in rows, samples in columns). 
#'   The input should use the same gene annotation system as the reference object.
#' @param xcell2object A pre-trained reference object of class \code{xCell2Object}, created using the \code{\link{xCell2Train}} function. 
#'   Pre-trained references are available within the package for common use cases.
#' @param minSharedGenes Minimum fraction of shared genes required between the mixture and the reference object (default: \code{0.9}). 
#'   If the shared fraction falls below this threshold, the function will stop with an error or warning, as accurate analysis 
#'   depends on sufficient overlap between the mixture and reference genes.
#' @param rawScores A Boolean indicating whether to return raw enrichment scores (default: \code{FALSE}). 
#'   Raw enrichment scores are computed directly from gene rankings without linear transformation or spillover correction.
#' @param spillover A Boolean to enable spillover correction on the enrichment scores (default: \code{TRUE}). 
#'   Spillover occurs when gene expression patterns overlap between closely related cell types, potentially inflating enrichment scores. 
#'   Correcting for spillover enhances the specificity of enrichment scores, particularly for related cell types. 
#'   The strength of this correction can be adjusted using the \code{spilloverAlpha} parameter.
#' @param spilloverAlpha A numeric value controlling the strength of spillover correction (default: \code{0.5}). 
#'   Lower values apply weaker correction, while higher values apply stronger correction. 
#'   An alpha value of 0.5 is suitable for most cases, but users may tune this parameter based on the similarity 
#'   of cell types in their reference.
#' @param BPPARAM A \linkS4class{BiocParallelParam} instance that determines the parallelization strategy (more in "Details"). 
#'   Default is \code{BiocParallel::SerialParam()}.
#' 
#' @return A data frame containing the cell type enrichment for each sample in the input matrix, as estimated by xCell2. 
#' Each row corresponds to a cell type, and each column corresponds to a sample.
#' 
#' @details
#' The \code{xCell2Analysis} function performs cell type enrichment analysis by leveraging gene signatures 
#' from a pre-trained \code{xCell2Object}. It computes enrichment scores for each cell type in the provided 
#' bulk gene expression mixture (\code{mix}), applies linear transformations, and corrects for spillover. 
#' Spillover correction addresses the overlap of gene expression patterns between closely related cell types, 
#' improving the specificity of the enrichment scores.
#'
#' ## Parallelization with \code{BPPARAM}
#' To achieve faster processing by running computations in parallel, \code{xCell2Analysis} supports parallelization through the \code{BPPARAM} 
#' parameter. Users can define a parallelization strategy using \code{BiocParallelParam} from the \code{BiocParallel} package. 
#' For example, \code{\link[BiocParallel]{MulticoreParam}} is suitable for multi-core processing on Linux and macOS, while 
#' \code{\link[BiocParallel]{SnowParam}} or \code{\link[BiocParallel]{SerialParam}} are better suited for Windows systems. 
#' Refer to the \href{https://www.bioconductor.org/packages/release/bioc/html/BiocParallel.html}{BiocParallel documentation} 
#' for further guidance on parallelization strategies.
#'
#' ## Relationship with Other Function(s)
#' The pre-trained \code{xCell2Object} used in \code{xCell2Analysis} is created via the \code{\link{xCell2Train}} function.
#' 
#' @examples
#' # For detailed example read xCell2 vignette.
#'
#' library(xCell2)
#'
#' # Load "ready to use" xCell2 reference object or generate a new one using `xCell2Train`
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#'
#' # Load demo bulk RNA-Seq gene expression mixture
#'
#' @param mix A bulk mixture of gene expression data (genes in rows, samples in columns). 
#'   The input must use the same gene annotation system as the reference object.
#' @param xcell2object A pre-trained reference object of class \code{xCell2Object}, created using \code{\link{xCell2Train}}. 
#'   Pre-trained references for common cases are provided within the package.
#' @param minSharedGenes Minimum fraction of shared genes required between the mixture and the reference object (default: \code{0.9}). 
#'   If the shared fraction is below this threshold, the function stops with an error or warning, as sufficient overlap is necessary 
#'   for accurate analysis.
#' @param rawScores Logical; if \code{TRUE}, returns raw enrichment scores without transformation or spillover correction (default: \code{FALSE}).
#' @param spillover Logical; enables spillover correction on enrichment scores (default: \code{TRUE}). 
#'   Spillover occurs when closely related cell types share gene expression patterns, inflating enrichment scores. 
#'   Correction enhances specificity, particularly for related cell types.
#' @param spilloverAlpha Numeric value controlling spillover correction strength (default: \code{0.5}). 
#'   Lower values apply weaker correction, while higher values apply stronger correction.
#' @param singscoreStableGenes description
#' @param BPPARAM A \linkS4class{BiocParallelParam} instance to define parallelization strategy (see "Details"). 
#'   Default is \code{BiocParallel::SerialParam()}.
#'
#' @return A data frame containing enrichment scores for each cell type and sample. 
#'   Rows correspond to cell types and columns to samples.
#'
#' @details
#' The \code{xCell2Analysis} function computes enrichment scores for each cell type using gene signatures 
#' from a pre-trained \code{xCell2Object}. Linear transformations and spillover corrections refine the results, 
#' improving specificity when cell types have overlapping gene expression patterns.
#'
#' \strong{Parallelization with \code{BPPARAM}:}
#' Computations can be parallelized using the \code{BPPARAM} parameter.  
#' Supported strategies include:
#' \itemize{
#'   \item \code{\link[BiocParallel]{MulticoreParam}} for multi-core processing (Linux/macOS).
#'   \item \code{\link[BiocParallel]{SnowParam}} or \code{\link[BiocParallel]{SerialParam}} for Windows systems.
#' }
#' See the \href{https://www.bioconductor.org/packages/release/bioc/html/BiocParallel.html}{BiocParallel documentation}.
#'
#' \strong{Relationship with Other Functions:}
#' The input reference object (\code{xCell2Object}) is created via \code{\link{xCell2Train}}.
#'
#' @examples
#' # For detailed examples, see the xCell2 vignette.
#'
#' library(xCell2)
#'
#' # Load pre-trained reference object
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#'
#' # Load demo bulk gene expression mixture
#' data(mix_demo, package = "xCell2")
#'
#' # Perform cell type enrichment analysis
#' xcell2_res <- xCell2::xCell2Analysis(
#'   mix = mix_demo, 
#'   xcell2object = DICE_demo.xCell2Ref
#' )
#'
#' # Parallel processing example with BiocParallel
#' library(BiocParallel)
#' parallel_param <- MulticoreParam(workers = 2)
#' xcell2_res_parallel <- xCell2::xCell2Analysis(
#'   mix = mix_demo, 
#'   xcell2object = DICE_demo.xCell2Ref, 
#'   BPPARAM = parallel_param
#' )
#'
#' # Example using stable genes for ranking in singscore
#' # Get a list of 30 stable genes (e.g., for cancer) from the singscore package
#' # Important: make sure those genes exist in your input mixture!
#' stable_genes <- singscore::getStableGenes(n_stable = 30, type = "carcinoma", id = "geneid")
#' 
#' # Note! Here we will use "fake" stable genes, just because the genes above not exist in "mix_demo"
#' fake_stable_genes <- rownames(mix_demo)[1:30]
#'
#' # Perform analysis using the stable genes for ranking
#' xcell2_res_stable <- xCell2::xCell2Analysis(
#'   mix = mix_demo,
#'   xcell2object = DICE_demo.xCell2Ref,
#'   singscoreStableGenes = fake_stable_genes
#' )
#' 
#' @seealso 
#' \code{\link{xCell2Train}}, for creating the reference object used in this analysis.
#'
#' @author Almog Angel and Dvir Aran
#' @export
xCell2Analysis <- function(mix,
                           xcell2object,
                           minSharedGenes = 0.9,
                           rawScores = FALSE,
                           spillover = TRUE,
                           spilloverAlpha = 0.5,
                           singscoreStableGenes = NULL,
                           BPPARAM = BiocParallel::SerialParam()) {
  
  message("Starting xCell2 Analysis...")
  
  scoreMix <- function(cellType, mixRanked, signaturesCellType) {

    scores <- vapply(signaturesCellType, function(sig) {
      singscore::simpleScore(mixRanked, upSet = sig, centerScore = FALSE)$TotalScore
    }, FUN.VALUE = double(ncol(mixRanked)))
    rownames(scores) <- colnames(mixRanked)
    return(scores)
  }
  calcEnrichment <- function(cellType) {
    
    signaturesCellType <- getSignatures(xcell2object)[startsWith(names(
      getSignatures(xcell2object)
    ), paste0(cellType, "#"))]
    
    # Remove signature with missing genes
    sigs2remove <- vapply(signaturesCellType, function(sig){
      sum(sig %in% rownames(mixRanked)) == 0
    }, FUN.VALUE = logical(1))
    
    if (all(sigs2remove) | sum(!sigs2remove) < 3) {
      warning(
        "Cannot calculate enrichment scores for '",
        cellType, "' because all signatgures' genes are missing in your mixture."
      )
      zeros_matrix <- matrix(rep(0, ncol(mixRanked)*3), nrow = ncol(mixRanked), ncol = 3)
      rownames(zeros_matrix) <- colnames(mixRanked)
      return(zeros_matrix)
    }
    
    if (any(sigs2remove)) {
      warning(
        "Removing ", sum(sigs2remove), "/", length(signaturesCellType), " signatures of '",
        cellType, "' because all genes in those signatures are missing in your mixture."
      )
    }
    
    signaturesCellType <- signaturesCellType[!sigs2remove]
    
    # Calculate enrichment scores with SingScore
    warning_list <- list()
    
    scores <- withCallingHandlers(
      {
        scoreMix(cellType, mixRanked, signaturesCellType)
      },
      warning = function(w) {
        warning_list <- c(warning_list, conditionMessage(w))  # Capture warnings
        invokeRestart("muffleWarning")  # Prevent the warning from being printed repeatedly
      }
    )
    
    if (length(warning_list) > 0) {
      warning("There were ",  length(warning_list), "/", length(signaturesCellType),
              " signatures with missing gene(s) for '", cellType,
              "' during enrichment scores calculation.\n",
              "Missing genes are:\n",
              unique(unlist(strsplit(gsub(".*genes missing: ", "", unlist(warning_list)), ","))))
    }
    
    return(scores)
  }
  
  # Wrapper function to include progress bar
  calcEnrichmentWrapper <- function(cellType) {
    pb$tick()
    calcEnrichment(cellType)
  }

  # Check reference/mixture genes intersection
  shared_genes <- intersect(rownames(mix), getGenesUsed(xcell2object))
  genesIntersectFrac <- round(length(shared_genes) / length(getGenesUsed(xcell2object)), 2)
  if (genesIntersectFrac < minSharedGenes) {
    stop(
      "This xCell2 reference shares ",
      genesIntersectFrac, " genes with the mixtures and minSharedGenes = ",
      minSharedGenes, ".",
      "\n", "Consider training a new xCell2 reference or adjusting minSharedGenes."
    )
  }

  if (genesIntersectFrac < 0.85) {
    warning(
      "This xCell2 reference shares only ",
      genesIntersectFrac, " genes with the mixtures.",
      "\n", "Consider using a reference that share at least 85% of the genes with the mixture for optimal results."
    )    
  }
  

  # Rank mix gene expression matrix
  if (is.null(singscoreStableGenes)) {
    mixRanked <- singscore::rankGenes(mix[shared_genes, ])
  }else{
    shared_stable_genes <- singscoreStableGenes[singscoreStableGenes %in% rownames(mix)]
    if (length(shared_stable_genes) == 0) {
      warning(
        "All stable genes not present in the mixture - Using singscore without stable genes normalization..."
      )
      mixRanked <- singscore::rankGenes(mix[shared_genes, ])
    }else{
      if (length(shared_stable_genes) < length(singscoreStableGenes)) {
        warning(
          "The following stable genes not present in the mixture: ",
          paste0(singscoreStableGenes[!singscoreStableGenes %in% rownames(mix)], collapse = " "),
          "\nUsing stable genes normalization with ",
          length(shared_stable_genes), " shared genes..."
        )
        mixRanked <- singscore::rankGenes(mix[shared_genes, ], stableGenes = shared_stable_genes)
      }else{
        message("Using stable genes to normalize ranks.")
        mixRanked <- singscore::rankGenes(mix[shared_genes, ], stableGenes = shared_stable_genes)
      }
    }
  }
  
  # Score and predict
  sigsCellTypes <- unique(unlist(lapply(
    names(getSignatures(xcell2object)),
    function(x) {
      strsplit(x, "#")[[1]][1]
    }
  )))
  
  # Get raw enrichment scores
  message("Calculating enrichment scores for all cell types...")
  
  if (BPPARAM$workers > 1) {
    resRaw <- BiocParallel::bplapply(sigsCellTypes, calcEnrichment, BPPARAM = BPPARAM)
  }else{
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent (:current/:total) :elapsed",
      total = length(sigsCellTypes), clear = FALSE, width = 60
    )
    resRaw <- BiocParallel::bplapply(sigsCellTypes, calcEnrichmentWrapper, BPPARAM = BPPARAM)
  }

  names(resRaw) <- sigsCellTypes
  
  res <- t(vapply(resRaw, function(cellTypeScores) {
    rowMeans(cellTypeScores)
  }, FUN.VALUE = double(nrow(resRaw[[1]]))))
  
  if (rawScores) {
    message("Returning raw enrichment scores without linear transformation or correction.")
    return(res)
  } else {
    # Linear transformation
    res <- t(vapply(rownames(res), function(cellType) {
      cellTypeRes <- res[cellType, ]
      
      # Linear transformation
      a <- getParams(xcell2object)[xcell2object@params$celltype == cellType, ]$a
      b <- getParams(xcell2object)[xcell2object@params$celltype == cellType, ]$b
      m <- getParams(xcell2object)[xcell2object@params$celltype == cellType, ]$m
      
      cellTypeRes <- (cellTypeRes^(1 / b)) / a
      cellTypeRes <- cellTypeRes * m
      
      # Shift values
      cellTypeRes <- cellTypeRes - min(cellTypeRes)
      cellTypeRes <- round(cellTypeRes, 5)
      
      return(cellTypeRes)
    }, FUN.VALUE = double(ncol(res))))
  }
  
  if (spillover) {
    message("Performing spillover correction...")
    # Spillover correction
    spillMat <- getSpillMat(xcell2object) * spilloverAlpha
    diag(spillMat) <- 1
    
    rows <- intersect(rownames(res), rownames(spillMat))
    
    scoresCorrected <- apply(res[rows, ], 2, function(x) {
      pracma::lsqlincon(spillMat[rows, rows],
                        x,
                        lb = 0
      )
    })
    scoresCorrected[scoresCorrected < 0] <- 0
    rownames(scoresCorrected) <- rows
    
    message("xCell2 Analysis completed successfully.")
    return(scoresCorrected)
  } else {
    message("Skipping spillover correction...")
    message("xCell2 Analysis completed successfully.")
    return(res)
  }
}
