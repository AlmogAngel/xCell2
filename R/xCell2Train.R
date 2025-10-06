# Function to validate input parameters for xCell2 training
ValidateInputs <- function(ref, labels, refType) {
  
  if (!refType %in% c("rnaseq", "array", "sc")) {
    stop("refType should be 'rnaseq', 'array' or 'sc'.")
  }
  
  # Check if input is a SummarizedExperiment or SingleCellExperiment
  if (inherits(ref, "SummarizedExperiment") || inherits(ref, "SingleCellExperiment")) {
    
    se <- ref
    
    if (is.null(labels)) {
      labels <- as.data.frame(colData(se))
    }else{
      if (!"data.frame" %in% class(labels)) {
        stop("labels must be a dataframe.")
      }
    }
    
    if (refType != "array") {
      # Attempt to access TPM data first, followed by other assay types
      if ("tpm" %in% assayNames(se)) {
        ref <- assay(se, "tpm")
      } else if ("logcounts" %in% assayNames(se)) {
        warning("TPM data not found. Using log-transformed counts instead.")
        ref <- assay(se, "logcounts")
      } else if ("normcounts" %in% assayNames(se)) {
        warning("TPM data not found. Using normalized counts instead.")
        ref <- assay(se, "normcounts")
      } else if ("counts" %in% assayNames(se)) {
        warning("TPM data not found. Using raw counts instead.")
        ref <- assay(se, "counts")
      } else {
        stop("No valid assay (tmp, logcounts, normcounts, counts) found in the SummarizedExperiment/SingleCellExperiment reference.")
      }
    } else {
      # Use counts for microarray data
      if ("counts" %in% assayNames(se)) {
        ref <- assay(se, "counts")
      } else {
        stop("No counts data found in the SummarizedExperiment object for microarray reference")
      }
    }
    
  } else {
    
    # Handle non-SummarizedExperiment/SingleCellExperiment input
    if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix", "dgRMatrix"))) {
      stop("ref must be one of these classes: matrix, dgCMatrix, Matrix, dgRMatrix or SummarizedExperiment/SingleCellExperiment object")
    }
    if (is.null(labels)) {
      stop("labels must be a data frame with 4 columns: 'ont', 'label', 'sample', and 'dataset'")
    }
    if (!"data.frame" %in% class(labels)) {
      stop("labels must be a dataframe.")
    }
    
  }
  
  if (all(c("ont", "label", "sample", "dataset") %in% colnames(labels))) {
    labels <- labels[, c("ont", "label", "sample", "dataset")]
  } else {
    if (inherits(ref, "SummarizedExperiment") || inherits(ref, "SingleCellExperiment")) {
      stop("colData(ref) must have 4 columns: 'ont'', 'label'', 'sample'' and 'dataset'")
      
    }else{
      stop("labels must have 4 columns: 'ont'', 'label'', 'sample'' and 'dataset'")
    }
  }
  
  if (length(unique(labels$label)) < 3) {
    stop("Reference must have at least 3 cell types!")
  }
  
  if (sum(grepl("_", labels$label)) != 0) {
    message("Changing underscores to dashes in cell-type labels!")
    labels$label <- gsub("_", "-", labels$label)
  }
  
  out <- list(ref = ref, labels = labels)
  return(out)
}

# Function to convert single-cell data to pseudo-bulk data
ScToPseudoBulk <- function(ref, labels, minPbCells, minPbSamples, BPPARAM) {
  cellTypes <- unique(labels$label)
  
  groupsList <- BiocParallel::bplapply(cellTypes, function(cellType) {
    cellTypeSamples <- labels[labels$label == cellType, ]$sample
    
    # Calculate maximum possible number of groups given minPbCells
    numGroups <- ceiling(length(cellTypeSamples) / minPbCells)
    
    if (length(cellTypeSamples) < minPbSamples) {
      minPbSamples <- length(cellTypeSamples)
    }
    
    if (numGroups < minPbSamples) {
      numGroups <- minPbSamples
    }
    
    # Generate minPbSamples pseudo samples of cellType
    if (length(cellTypeSamples) > minPbSamples) {
      cellTypeSamplesShuffled <- sample(cellTypeSamples, length(cellTypeSamples))
      listOfShuffledSamples <- split(
        cellTypeSamplesShuffled,
        ceiling(seq_along(cellTypeSamplesShuffled)
                / (length(cellTypeSamplesShuffled) / numGroups))
      )
      
      tmp <- vapply(listOfShuffledSamples, function(group) {
        if (length(group) == 1) {
          ref[, group]
        } else {
          if ("matrix" %in% class(ref)) {
            Rfast::rowsums(ref[, group])
          } else {
            Matrix::rowSums(ref[, group])
          }
        }
      }, FUN.VALUE = double(nrow(ref)))
    } else {
      tmp <- ref[, cellTypeSamples]
      tmp <- as.data.frame(as.matrix(tmp))
    }
    
    colnames(tmp) <- paste0(cellType, ".", seq_len(ncol(tmp)))
    tmp
    
  }, BPPARAM = BPPARAM)
  names(groupsList) <- cellTypes
  
  pseudoRef <- as.matrix(Reduce(cbind, groupsList))
  rownames(pseudoRef) <- rownames(ref)
  
  label_tmp <- tibble::tibble(
    label = sub("\\.\\d+$", "", colnames(pseudoRef)),
    sample = colnames(pseudoRef),
    dataset = "pseudoBulk"
  )
  
  pseudoLabel <- label_tmp %>%
    left_join(unique(labels[,c("ont", "label")]), by = "label") %>%
    dplyr::select(ont, everything()) %>%
    as.data.frame()
  
  return(list(ref = pseudoRef, labels = pseudoLabel))
}

# Function to prepare reference and mixture data
PrepRefMix <- function(ref, mix, refType, minScGenes) {
  
  if (refType == "sc") {
    message("> Normalizing pseudo bulk reference to CPM.")
    libSizes <- Matrix::colSums(ref)
    normFactor <- 1e6 / libSizes
    refNorm <- ref %*% Matrix::Diagonal(x = normFactor)
    colnames(refNorm) <- colnames(ref)
    ref <- as.matrix(refNorm)
    message("> Filtering pseudo bulk genes by variance.")
    genesVar <- sort(apply(ref, 1, var), decreasing = TRUE)
    varCutoff <- c(1.5, 1, 0.8, 0.5, 0.3, 0.1, 0)
    for (co in varCutoff) {
      if (co == 0) {
        varGenesToUse <- names(genesVar)[seq_len(round(length(genesVar) * 0.5))]
      } else {
        varGenesToUse <- names(genesVar[genesVar >= co])
      }
      if (length(varGenesToUse) > minScGenes) {
        break
      } else {
        stop("Not enough variable genes in scRNA-Seq reference!")
      }
    }
    ref <- ref[varGenesToUse, ]
    
    if (min(ref) < 3) {
      ref <- ref + 3
    }
  } else {
    if (max(ref) < 50) {
      ref <- 2^ref
      if (min(ref) == 1) {
        ref <- ref - 1
      }
      
      if (min(ref) < 3) {
        ref <- ref + 3
      }
    } else {
      if (min(ref) < 3) {
        ref <- ref + 3
      }
    }
  }
  
  ref <- log2(ref)
  
  if (!is.null(mix)) {
    if (max(mix) >= 50) {
      mix <- log2(mix + 1)
    }
    
    sharedGenes <- intersect(rownames(ref), rownames(mix))
    wm <- paste0(length(sharedGenes), " genes are shared between reference and mixture.")
    message(sprintf(wm))
    ref <- ref[sharedGenes, ]
    mix <- mix[sharedGenes, ]
  }
  
  return(list(refOut = ref, mixOut = mix))
}

# Function to get cell type dependencies from a lineage file
LoadDependenciesFromFile <- function(lineageFileChecked) {
  ont <- readr::read_tsv(lineageFileChecked, show_col_types = FALSE) %>%
    dplyr::mutate_all(as.character)
  
  cellTypes <- dplyr::pull(ont[, 2])
  cellTypes <- gsub("_", "-", cellTypes)
  depList <- vector(mode = "list", length = length(cellTypes))
  names(depList) <- cellTypes
  
  for (i in seq_len(nrow(ont))) {
    descendants <- gsub("_", "-", strsplit(dplyr::pull(ont[i, 3]), ";")[[1]])
    descendants <- descendants[!is.na(descendants)]
    
    ancestors <- gsub("_", "-", strsplit(dplyr::pull(ont[i, 4]), ";")[[1]])
    ancestors <- ancestors[!is.na(ancestors)]
    
    depList[[i]] <- list("descendants" = descendants, "ancestors" = ancestors)
  }
  
  return(depList)
}

# Function to calculate quantiles for each cell type
MakeQuantiles <- function(ref, labels, probs, BPPARAM) {
  
  cellTypes <- unique(labels[, 2])
  
  quantilesMatList <- BiocParallel::bplapply(cellTypes, function(type) {
    typeSamples <- labels[, 2] == type
    
    if (sum(typeSamples) == 1) {
      type.df <- cbind(ref[, typeSamples], ref[, typeSamples])
    } else {
      type.df <- ref[, typeSamples]
    }
    
    # Calculate quantiles
    apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1 - probs))), na.rm = TRUE))
  }, BPPARAM = BPPARAM)
  names(quantilesMatList) <- cellTypes
  
  return(quantilesMatList)
}

# Function to create signatures for each cell type
CreateSignatures <- function(labels, depList, quantilesMatrix, probs, diffVals, minGenes,
                             maxGenes, minFracCtPassed, BPPARAM) {
  # Inner function to generate signatures for a single cell type
  getSignatures <- function(cellTypes, type, depList, quantilesMatrix, probs,
                            diffVals, minGenes, maxGenes, minFracCtPassed) {
    # Remove dependent cell types
    notDepCellTypes <- cellTypes[!cellTypes %in% c(type, unname(unlist(depList[[type]])))]
    
    # Set signature cutoffs grid
    paramDf <- expand.grid("diffVals" = diffVals, "probs" = probs)
    
    # Find top genes
    typeSignatures <- list()
    
    while (length(typeSignatures) < 3 & minFracCtPassed >= 0) {
      maxGenesProblem <- c()
      for (i in seq_len(nrow(paramDf))) {
        # Get Boolean matrices with genes that pass the quantiles criteria
        diff <- paramDf[i, ]$diffVals
        lowerProb <- which(as.character(as.numeric(gsub(
          "%", "", rownames(quantilesMatrix[[1]])
        )) / 100) == as.character(paramDf[i, ]$probs))
        upperProb <- which(as.character(as.numeric(gsub(
          "%", "", rownames(quantilesMatrix[[1]])
        )) / 100) == as.character(1 - paramDf[i, ]$probs))
        
        upperProbMatrix <- vapply(notDepCellTypes, function(x) {
          get(x, quantilesMatrix)[upperProb, ]
        }, FUN.VALUE = double(ncol(quantilesMatrix[[1]])))
        
        # Check diff-prob criteria
        diffGenesMatrix <- apply(upperProbMatrix, 2, function(x) {
          get(type, quantilesMatrix)[lowerProb, ] > x + diff
        })
        
        genesScores <- apply(diffGenesMatrix, 1, function(x) {
          names(which(x))
        })
        genesScores <- genesScores[lengths(genesScores) > 0]
        nCtPassed <- sort(unique(lengths(genesScores)), decreasing = TRUE)
        
        # Make signatures
        for (j in nCtPassed) {
          fracCtPassed <- round(j / length(notDepCellTypes), 2)
          if (fracCtPassed < minFracCtPassed) break
          
          sigGenes <- names(which(lengths(genesScores) >= j))
          nGenes <- length(sigGenes)
          
          if (nGenes < minGenes) next
          if (nGenes > maxGenes) {
            maxGenesProblem <- c(maxGenesProblem, TRUE)
            break
          } else {
            maxGenesProblem <- c(maxGenesProblem, FALSE)
          }
          
          sigName <- paste(paste0(type, "#"), paramDf[i, ]$probs, diff,
                           nGenes, fracCtPassed,
                           sep = "_"
          )
          typeSignatures[[sigName]] <- sigGenes
        }
      }
      
      # Handle cases where too many genes are differentially expressed
      if (all(maxGenesProblem)) {
        diffValsStrict <- c(
          diffVals, log2(2^max(diffVals) * 2),
          log2(2^max(diffVals) * 4),
          log2(2^max(diffVals) * 8),
          log2(2^max(diffVals) * 16),
          log2(2^max(diffVals) * 32)
        )
        paramDf <- expand.grid("diffVals" = diffValsStrict, "probs" = probs)
        
        for (i in seq_len(nrow(paramDf))) {
          lowerProb <- which(as.character(as.numeric(gsub(
            "%", "", rownames(quantilesMatrix[[1]])
          )) / 100) == as.character(paramDf[i, ]$probs))
          upperProb <- which(as.character(as.numeric(gsub(
            "%", "", rownames(quantilesMatrix[[1]])
          )) / 100) == as.character(
            1 - paramDf[i, ]$probs
          ))
          
          upperProbMatrix <- vapply(notDepCellTypes, function(x) {
            get(x, quantilesMatrix)[upperProb, ]
          }, FUN.VALUE = double(ncol(quantilesMatrix[[1]])))
          
          
          diffGenesMatrix <- apply(upperProbMatrix, 2, function(x) {
            get(type, quantilesMatrix)[lowerProb, ] > x + diff
          })
          
          genesScores <- apply(diffGenesMatrix, 1, function(x) {
            names(which(x))
          })
          genesScores <- genesScores[lengths(genesScores) > 0]
          nCtPassed <- sort(unique(lengths(genesScores)), decreasing = TRUE)
          
          # Make signatures
          for (j in nCtPassed) {
            fracCtPassed <- round(j / length(notDepCellTypes), 2)
            if (fracCtPassed < minFracCtPassed) break
            
            sigGenes <- names(which(lengths(genesScores) >= j))
            nGenes <- length(sigGenes)
            
            if (nGenes < minGenes) next
            if (nGenes > maxGenes) {
              maxGenesProblem <- c(maxGenesProblem, TRUE)
              break
            }
            
            sigName <- paste(paste0(type, "#"), paramDf[i, ]$probs,
                             diff, nGenes, fracCtPassed,
                             sep = "_"
            )
            typeSignatures[[sigName]] <- sigGenes
          }
        }
      }
      
      # Remove duplicate signatures
      typeSignaturesSorted <- lapply(typeSignatures, function(x) sort(x))
      typeSignaturesCollapsed <- vapply(typeSignaturesSorted,
                                        paste,
                                        collapse = ",", FUN.VALUE = character(1)
      )
      duplicatedSignatures <- duplicated(typeSignaturesCollapsed)
      typeSignatures <- typeSignatures[!duplicatedSignatures]
      
      # Relax parameter until there are at least 3 signatures
      minFracCtPassed <- minFracCtPassed - 0.05
    }
    
    return(typeSignatures)
  }
  
  
  cellTypes <- unique(labels[, 2])
  
  allSignatures <- BiocParallel::bplapply(cellTypes, function(type) {
    typeSignatures <- getSignatures(
      cellTypes, type, depList, quantilesMatrix, probs,
      diffVals, minGenes, maxGenes, minFracCtPassed
    )
    
    # Ensure minimum 3 signatures per cell type
    if (length(typeSignatures) < 3) {
      probsToUse <- c(probs, max(probs) * 1.25, max(probs) * 1.5, max(probs) * 2)
      probsToUse <- probsToUse[probsToUse < 0.5]
      minDiff <- min(diffVals[diffVals != 0])
      diffValsToUse <- unique(c(minDiff * 0, minDiff * 0.5, minDiff * 0.75, diffVals))
      minGenesToUse <- round(minGenes * 0.5)
      minGenesToUse <- ifelse(minGenesToUse < 3, 3, minGenesToUse)
      
      typeSignatures <- getSignatures(cellTypes, type, depList, quantilesMatrix,
                                      probs = probsToUse, diffVals = diffValsToUse,
                                      minGenes = minGenesToUse, maxGenes, minFracCtPassed
      )
    }
    
    return(typeSignatures)
  }, BPPARAM = BPPARAM)
  
  allSignatures <- unlist(allSignatures, recursive = FALSE)
  
  if (length(allSignatures) == 0) {
    stop("No signatures found for reference!")
  }
  
  return(allSignatures)
}

# Function to generate the gene expression profile matrix
MakeGEPMat <- function(ref, labels) {
  cellTypes <- unique(labels$label)
  
  gepMat <- vapply(cellTypes, function(type) {
    typeSamples <- labels[, 2] == type
    if (sum(typeSamples) == 1) {
      typeVec <- as.vector(ref[, typeSamples])
    } else {
      typeVec <- Rfast::rowMedians(as.matrix(ref[, typeSamples]))
    }
  }, FUN.VALUE = double(nrow(ref)))
  rownames(gepMat) <- rownames(ref)
  
  return(gepMat)
}

# Function to calculate cell types correlation
GetCellTypeCorrelation <- function(gepMat, refType) {
  cellTypes <- colnames(gepMat)
  
  if (refType != "sc") {
    # Use top 10% most variable genes
    genesVar <- apply(gepMat, 1, var)
    mostVarGenesCutoff <- quantile(genesVar, 0.9, na.rm = TRUE)
    gepMat <- gepMat[genesVar > mostVarGenesCutoff, ]
  } else {
    # Use top 1% most variable genes
    genesVar <- apply(gepMat, 1, var)
    mostVarGenesCutoff <- quantile(genesVar, 0.99, na.rm = TRUE)
    gepMat <- gepMat[genesVar > mostVarGenesCutoff, ]
  }
  
  # Create correlation matrix
  corMat <- matrix(1, ncol = length(cellTypes), nrow = length(cellTypes), dimnames = list(cellTypes, cellTypes))
  lowerTriCoord <- which(lower.tri(corMat), arr.ind = TRUE)
  
  # Calculate correlations
  for (i in seq_len(nrow(lowerTriCoord))) {
    celltypeI <- rownames(corMat)[lowerTriCoord[i, 1]]
    celltypeJ <- colnames(corMat)[lowerTriCoord[i, 2]]
    corMat[lowerTriCoord[i, 1], lowerTriCoord[i, 2]] <- cor(gepMat[, celltypeI], gepMat[, celltypeJ], method = "spearman")
    corMat[lowerTriCoord[i, 2], lowerTriCoord[i, 1]] <- cor(gepMat[, celltypeI], gepMat[, celltypeJ], method = "spearman")
  }
  
  return(corMat)
}

# Function to learn linear transformation and spillover parameters
LearnParams <- function(gepMat, corMat, signatures, depList, BPPARAM) {
  
  
  cellTypes <- colnames(gepMat)
  gepMatLinear <- 2^gepMat
  simFracs <- c(0, seq(0.01, 0.25, 0.01))
  fracToUse <- 0.25
  topSpillValue <- 0.5
  
  # Generate mixtures
  mixList <- BiocParallel::bplapply(cellTypes, function(cellType) {
    # Generate cellType mixture
    cellTypeMat <- matrix(rep(gepMatLinear[, cellType], length(simFracs)),
                          byrow = FALSE, ncol = length(simFracs),
                          dimnames = list(rownames(gepMatLinear), simFracs)
    )
    cellTypeMatFrac <- cellTypeMat %*% diag(simFracs)
    
    # Generate control mixture
    if (!is.null(depList)) {
      depCts <- unique(c(cellType, unname(unlist(depList[[cellType]]))))
      controls <- cellTypes[!cellTypes %in% depCts]
    } else {
      controls <- cellTypes[cellTypes != cellType]
    }
    
    control <- names(sort(corMat[controls, cellType])[1])
    controlsMat <- matrix(rep(gepMatLinear[, control], length(simFracs)),
                          byrow = FALSE, ncol = length(simFracs),
                          dimnames = list(rownames(gepMatLinear), simFracs)
    )
    controlsMatFrac <- controlsMat %*% diag(1 - simFracs)
    
    # Combine
    mixture <- cellTypeMatFrac + controlsMatFrac
    colnames(mixture) <- paste0(cellType, "^^", control, "%%", simFracs)
    
    return(mixture)
  }, BPPARAM = BPPARAM)
  names(mixList) <- cellTypes
  
  # Learn linear parameters
  linearParams <- BiocParallel::bplapply(cellTypes, function(cellType) {
    # Get scores
    mixMatRanked <- singscore::rankGenes(mixList[[cellType]])
    signaturesCellType <- signatures[gsub("#.*", "", names(signatures)) %in% cellType]
    scores <- rowMeans(vapply(signaturesCellType, function(sig) {
      singscore::simpleScore(mixMatRanked, upSet = sig, centerScore = FALSE)$TotalScore
    }, FUN.VALUE = double(ncol(mixMatRanked))))
    
    # Get transformation parameters
    tp <- try(minpack.lm::nlsLM(scores ~ a * simFracs^b,
                                start = list(a = 1, b = 1),
                                control = list(maxiter = 500)
    ), silent = TRUE)
    a <- coef(tp)[[1]]
    b <- coef(tp)[[2]]
    
    # Get linear model parameters
    scoresTransformed <- (scores^(1 / b)) / a
    lmFit <- lm(simFracs ~ scoresTransformed)
    m <- coef(lmFit)[[2]]
    n <- coef(lmFit)[[1]]
    
    return(tibble::tibble("celltype" = cellType, "a" = a, "b" = b, "m" = m, "n" = n))
  }, BPPARAM = BPPARAM)
  linearParams <- dplyr::bind_rows(linearParams)
  
  # Learn spillover parameters
  spillScores <- BiocParallel::bplapply(cellTypes, function(cellType) {
    signaturesCellType <- signatures[gsub("#.*", "", names(signatures)) %in% cellType]
    
    a <- linearParams[linearParams$celltype == cellType, ]$a
    b <- linearParams[linearParams$celltype == cellType, ]$b
    m <- linearParams[linearParams$celltype == cellType, ]$m
    n <- linearParams[linearParams$celltype == cellType, ]$n
    
    # Generate fracToUse mixtures
    ctsMatFrac <- gepMatLinear * fracToUse
    
    # Generate control mixture
    if (!is.null(depList)) {
      depCts <- unique(c(cellType, unname(unlist(depList[[cellType]]))))
      controls <- cellTypes[!cellTypes %in% depCts]
    } else {
      controls <- cellTypes[cellTypes != cellType]
    }
    
    controls <- vapply(colnames(ctsMatFrac), function(ct) {
      names(sort(corMat[controls, ct])[1])
    }, FUN.VALUE = character(1))
    controlsMatFrac <- gepMatLinear[, controls] * (1 - fracToUse)
    
    # Combine
    mixture <- ctsMatFrac + controlsMatFrac
    
    # Get results for all cell type mixtures
    mixCtsMatRanked <- singscore::rankGenes(mixture)
    mixCtsMatScores <- vapply(signaturesCellType, function(sig) {
      singscore::simpleScore(mixCtsMatRanked, upSet = sig, centerScore = FALSE)$TotalScore
    }, FUN.VALUE = double(ncol(mixCtsMatRanked)))
    mixCtsMatScores <- Rfast::rowmeans(mixCtsMatScores)
    mixCtsMatScores <- (mixCtsMatScores^(1 / b)) / a
    mixCtsMatScores <- mixCtsMatScores * m + n
    names(mixCtsMatScores) <- colnames(mixCtsMatRanked)
    
    # Get results for all cell type controls
    controlsCtsMatRanked <- singscore::rankGenes(controlsMatFrac)
    colnames(controlsCtsMatRanked) <- make.unique(colnames(controlsCtsMatRanked))
    controlsCtsMatScores <- vapply(signaturesCellType, function(sig) {
      singscore::simpleScore(controlsCtsMatRanked, upSet = sig, centerScore = FALSE)$TotalScore
    }, FUN.VALUE = double(ncol(controlsCtsMatRanked)))
    controlsCtsMatScores[controlsCtsMatScores < 0] <- 0
    controlsCtsMatScores <- Rfast::rowmeans(controlsCtsMatScores)
    controlsCtsMatScores <- (controlsCtsMatScores^(1 / b)) / a
    controlsCtsMatScores <- controlsCtsMatScores * m + n
    controlsCtsMatScores[controlsCtsMatScores < 0] <- 0
    names(controlsCtsMatScores) <- controls
    
    finalScores <- round(mixCtsMatScores - controlsCtsMatScores, 2)
    if (!is.null(depList)) {
      depCts <- depCts[depCts != cellType]
      finalScores[depCts] <- 0
    }
    finalScores[finalScores < 0] <- 0
    
    return(finalScores)
  }, BPPARAM = BPPARAM)
  names(spillScores) <- cellTypes
  
  # Create spillover matrix
  spillMat <- Reduce(rbind, spillScores)
  rownames(spillMat) <- cellTypes
  
  # Normalize and clean spillover matrix
  spillMat[spillMat < 0] <- 0
  spillMat <- spillMat / diag(spillMat)
  
  spillMat[is.nan(spillMat)] <- 0
  spillMat[spillMat > 1] <- 1
  
  spillMat[spillMat > topSpillValue] <- topSpillValue
  diag(spillMat) <- 1
  
  return(list(params = linearParams, spillmat = spillMat))
}


#' Train Custom xCell2 Reference Object
#'
#' This function creates a custom reference object for \code{\link{xCell2Analysis}}, enabling cell type enrichment analysis.
#' It supports references derived from RNA-Seq, microarray, and scRNA-Seq data and can be derived from various tissues and organisms. 
#'
#' @importFrom magrittr %>%
#' @importFrom tidyselect everything
#' @importFrom dplyr bind_cols select right_join mutate_all pull bind_rows left_join select
#' @importFrom tibble tibble
#' @importFrom readr read_tsv
#' @importFrom BiocParallel SerialParam MulticoreParam bplapply
#' @importFrom utils data
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SingleCellExperiment colData
#' @importFrom minpack.lm nlsLM
#' @importFrom Rfast rowMedians rowmeans rowsums Sort
#' @importFrom Matrix rowMeans rowSums colSums Diagonal
#' @importFrom singscore rankGenes simpleScore
#' @importFrom stats coef cor lm quantile var
#' @importFrom methods new
#'
#' @param ref A reference gene expression matrix (genes in rows, samples/cells in columns) or a 
#'   \linkS4class{SummarizedExperiment}/\linkS4class{SingleCellExperiment} object with expression data in the assays slot.  
#'   
#'   \strong{Valid Assays:}
#'   \describe{
#'     \item{\code{"tpm"}}{Transcripts Per Million (recommended for RNA-Seq).}
#'     \item{\code{"logcounts"}}{Log-transformed normalized counts.}
#'     \item{\code{"normcounts"}}{Normalized counts.}
#'     \item{\code{"counts"}}{Raw counts (required for microarray references).}
#'   }
#'   
#'   \strong{Notes:}
#'   \itemize{
#'     \item If is RNA-Seq data - normalization by gene length is highly recommended.
#'     \item If multiple assays exist, \code{"tpm"} is prioritized.
#'     \item For microarray data, the \code{"counts"} assay must be used.
#'   }
#'
#' @param mix A bulk mixture of gene expression matrix (genes in rows, samples in columns) (optional). 
#'   This parameter is required if \code{returnAnalysis} is set to \code{TRUE}, as it is used for enrichment analysis.
#'
#' @param labels A data frame with the following columns:
#'   \itemize{
#'     \item \code{"ont"}: The cell type ontology ID (e.g., \code{"CL:0000545"}). Set to \code{NA} if not available. 
#'       Ontologies can be found at \href{https://www.ebi.ac.uk/ols4/ontologies/cl}{EBI Ontology Lookup Service (OLS)} or 
#'       by using the \link[ontologyIndex]{ontologyIndex} package.
#'     \item \code{"label"}: The cell type name (e.g., \code{"T-helper 1 cell"}).
#'     \item \code{"sample"}: The sample or cell identifier, matching column names in the reference matrix.
#'     \item \code{"dataset"}: The dataset source for each sample. If not applicable, use a constant value for all samples.
#'   }
#'   If \code{ref} is a \code{SummarizedExperiment} or \code{SingleCellExperiment} object, this parameter can be used to 
#'   **override** the default labels extracted from \code{colData(ref)}.
#'   
#' @param refType The type of reference data: \code{"rnaseq"} for RNA-Seq, \code{"array"} for microarray, or \code{"sc"} for scRNA-Seq.
#'
#' @param minPbCells Minimum number of cells in a pseudo-bulk sample for scRNA-Seq references (default: \code{30}).
#' @param minPbSamples Minimum number of pseudo-bulk samples for scRNA-Seq references (default: \code{10}).
#' @param minScGenes Minimum number of genes for pseudo-bulk samples for scRNA-Seq references (default: \code{1e4}).
#'
#' @param useOntology A Boolean indicating whether to use ontological integration for cell type dependencies (default: \code{TRUE}). 
#'   Lineage relationships are determined using the Cell Ontology (CL). Users can refine these dependencies with 
#'   \code{\link{xCell2GetLineage}} and provide them via the \code{lineageFile} parameter.
#'
#' @param lineageFile Path to a manually curated cell type lineage file generated with \code{\link{xCell2GetLineage}} (optional).
#'
#' @param BPPARAM A \linkS4class{BiocParallelParam} instance that determines the parallelization strategy (more in "Details"). 
#'   Default is \code{BiocParallel::SerialParam()}.
#'
#' @param returnSignatures A Boolean to return only cell type signatures (default: \code{FALSE}).
#' @param returnAnalysis A Boolean to return \code{\link{xCell2Analysis}} results instead of a reference object (default: \code{FALSE}).
#'
#' @param useSpillover A Boolean to use spillover correction during analysis when \code{returnAnalysis} is \code{TRUE} (default: \code{TRUE}).
#'   Spillover correction enhances the specificity of enrichment scores by accounting for overlaps between cell types.
#'
#' @param spilloverAlpha Numeric value controlling spillover correction strength (default: \code{0.5}).
#'   Lower values apply weaker correction, while higher values apply stronger correction.
#'
#' @return An \code{xCell2Object} containing:
#'   \itemize{
#'     \item \strong{signatures}: Cell type-specific gene signatures.
#'     \item \strong{dependencies}: Lineage-based dependencies.
#'     \item \strong{params}: Linear transformation parameters.
#'     \item \strong{spill_mat}: A spillover correction matrix.
#'     \item \strong{genes_used}: Genes used for training.
#'   }
#'
#' @details
#' \strong{Ontological Integration:}
#' Ontological integration (\code{useOntology}) leverages hierarchical cell type relationships to ensure biologically meaningful signatures. 
#' Dependencies can be refined using \code{\link{xCell2GetLineage}}, which generates lineage files for manual review.
#'
#' \strong{Spillover Correction:}
#' Spillover correction enhances the specificity of enrichment scores by reducing overlaps between related cell types. 
#' Use the \code{spilloverAlpha} parameter to tune the strength of correction.
#'
#' \strong{Contribute Your xCell2 Reference Object:}
#' Users are encouraged to share their reference objects via the \href{https://dviraran.github.io/xCell2ref}{xCell2 Reference Repository}.
#'
#' @examples
#' library(xCell2)
#' data(dice_demo_ref, package = "xCell2")
#' dice_ref <- SummarizedExperiment::assay(dice_demo_ref, "logcounts")
#' colnames(dice_ref) <- make.unique(colnames(dice_ref))
#' dice_labels <- as.data.frame(SummarizedExperiment::colData(dice_demo_ref))
#' dice_labels$ont <- NA
#' dice_labels$sample <- colnames(dice_ref)
#' dice_labels$dataset <- "DICE"
#' DICE.xCell2Ref <- xCell2::xCell2Train(ref = dice_ref, labels = dice_labels, refType = "rnaseq")
#'
#' # Parallel processing example with BiocParallel
#' library(BiocParallel)
#' parallel_param <- MulticoreParam(workers = 2)
#' DICE.xCell2Ref <- xCell2::xCell2Train(ref = dice_ref, labels = dice_labels, refType = "rnaseq",
#'  BPPARAM = parallel_param)
#' 
#' @seealso 
#' \code{\link{xCell2Analysis}}, for enrichment analysis. 
#' \code{\link{xCell2GetLineage}}, for refining cell type dependencies.
#' @author Almog Angel and Dvir Aran
#' @export
xCell2Train <- function(ref,
                        mix = NULL,
                        labels = NULL,
                        refType,
                        lineageFile = NULL,
                        BPPARAM = BiocParallel::SerialParam(),
                        useOntology = TRUE,
                        returnSignatures = FALSE,
                        returnAnalysis = FALSE,
                        useSpillover = TRUE,
                        spilloverAlpha = 0.5,
                        minPbCells = 30,
                        minPbSamples = 10,
                        minScGenes = 1e4
) {
  
  message("Starting xCell2 Train...")
  
  inputsValidated <- ValidateInputs(ref, labels, refType)
  ref <- inputsValidated$ref
  labels <- inputsValidated$labels
  
  # Generate pseudo bulk from scRNA-Seq reference
  if (refType == "sc") {
    message("Converting scRNA-seq reference to pseudo bulk...")
    pbData <- ScToPseudoBulk(ref, labels, minPbCells, minPbSamples, BPPARAM)
    ref <- pbData$ref
    labels <- pbData$labels
  }
  
  # Prepare reference and mixture data
  out <- PrepRefMix(ref, mix, refType, minScGenes)
  ref <- out$refOut
  mix <- out$mixOut
  sharedGenes <- rownames(ref)
  
  # Get cell type dependencies list
  if (useOntology) {
    if (is.null(lineageFile)) {
      message("Finding dependencies using cell type ontology...")
      depList <- xCell2::xCell2GetLineage(labels, outFile = NULL)
    } else {
      message("Loading cell type dependencies...")
      depList <- LoadDependenciesFromFile(lineageFile)
    }
  } else {
    message("Skipping ontological integration!")
    depList <- NULL
  }
  
  # Generate signatures
  probs <- c(0.1, 0.25, 0.333, 0.49)
  diffVals <- c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5)
  minGenes <- 8
  maxGenes <- 200
  minFracCtPassed <- 0.5
  message("Generating signatures...")
  quantilesMatrix <- MakeQuantiles(ref, labels, probs, BPPARAM)
  signatures <- CreateSignatures(
    labels, depList, quantilesMatrix, probs, diffVals, minGenes,
    maxGenes, minFracCtPassed, BPPARAM
  )
  
  if (returnSignatures) {
    message("Retuning xCell2 reference object with signatures only.")
    xCell2S4 <- new("xCell2Object",
                    signatures = signatures,
                    dependencies = list(),
                    params = data.frame(),
                    spill_mat = matrix(),
                    genes_used = sharedGenes
    )
    return(xCell2S4)
  }
  
  # Learn linear transformation parameters
  message("Learning linear transformation and spillover parameters...")
  gepMat <- MakeGEPMat(ref, labels)
  corMat <- GetCellTypeCorrelation(gepMat, refType)
  params <- LearnParams(gepMat, corMat, signatures, depList, BPPARAM)
  
  # Save results in S4 object
  if (is.null(depList)) {
    depList <- list()
  }
  xCell2S4 <- new("xCell2Object",
                  signatures = signatures,
                  dependencies = depList,
                  params = params$params,
                  spill_mat = params$spillmat,
                  genes_used = sharedGenes
  )
  
  message("Your custom xCell2 reference object is ready!")
  message("> Please consider sharing with others here: https://dviraran.github.io/xCell2ref")
  
  if (returnAnalysis) {
    message("Running xCell2Analysis...")
    res <- xCell2::xCell2Analysis(mix,
                                  xcell2object = xCell2S4, spillover = useSpillover,
                                  spilloverAlpha = spilloverAlpha, BPPARAM = BPPARAM
    )
    return(res)
  } else {
    return(xCell2S4)
  }
}
