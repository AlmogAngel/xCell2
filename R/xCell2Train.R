# Function to validate input parameters for xCell2 training
ValidateInputs <- function(ref, labels, refType) {
  
  # Check if input is a SummarizedExperiment or SingleCellExperiment
  if (inherits(ref, "SummarizedExperiment") || inherits(ref, "SingleCellExperiment")) {
    se <- ref
    ref <- assays(se)$counts
    labels <- as.data.frame(colData(se))
  }else{
    if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
      stop("ref must be one of these classes: matrix, dgCMatrix, Matrix or SummarizedExperiment/SingleCellExperiment object")
    }
    if (is.null(labels)) {
      stop("labels must be a data frame with 4 columns: 'ont'', 'label'', 'sample'' and 'dataset'")
    }
  }
  
  if (all(colnames(labels) %in% c("ont", "label", "sample", "dataset"))) {
    labels <- labels[, c("ont", "label", "sample", "dataset")]
  } else {
    stop("labels must have 4 columns: 'ont'', 'label'', 'sample'' and 'dataset'")
  }
  
  if (length(unique(labels$label)) < 3) {
    stop("Reference must have at least 3 cell types!")
  }
  
  if (!"data.frame" %in% class(labels)) {
    stop("labels must be a dataframe.")
  }
  
  if (!refType %in% c("rnaseq", "array", "sc")) {
    stop("refType should be 'rnaseq', 'array' or 'sc'.")
  }
  
  if (sum(grepl("_", labels$label)) != 0) {
    message("Changing underscores to dashes in cell-type labels!")
    labels$label <- gsub("_", "-", labels$label)
  }
  
  out <- list(ref = ref, labels = labels)
  return(out)
}

# Function to convert single-cell data to pseudo-bulk data
ScToPseudoBulk <- function(ref, labels, minPbCells, minPbSamples) {
  cellTypes <- unique(labels$label)
  
  groupsList <- lapply(cellTypes, function(cellType) {
    cellTypeSamples <- labels[labels$label == cellType, ]$sample
    
    # Calculate maximum possible number of groups given minPbCells
    numGroups <- ceiling(length(cellTypeSamples) / minPbCells)
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
      
      vapply(listOfShuffledSamples, function(group) {
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
      colnames(tmp) <- as.character(seq_len(tmp))
      tmp
    }
  })
  names(groupsList) <- cellTypes
  
  pseudoRef <- as.matrix(dplyr::bind_cols(groupsList))
  rownames(pseudoRef) <- rownames(ref)
  
  pseudoLabel <- tibble::tibble(labels) %>%
    dplyr::select("ont", "label") %>%
    unique() %>%
    dplyr::right_join(tibble::tibble(
      label = sub("\\.\\d+$", "", colnames(pseudoRef)),
      sample = colnames(pseudoRef),
      dataset = "pseudoBulk"
    ), by = "label") %>%
    as.data.frame()
  
  return(list(ref = pseudoRef, labels = pseudoLabel))
}

# Function to prepare reference and mixture data
PrepRefMix <- function(ref, mix, refType, minScGenes, humanToMouse) {
  if (humanToMouse) {
    message("Converting reference genes from human to mouse...")
    local_env <- new.env()
    data(human_mouse_gene_symbols, package = "xCell2", envir = local_env)
    human_mouse_gene_symbols <- local_env$human_mouse_gene_symbols
    rownames(human_mouse_gene_symbols) <- human_mouse_gene_symbols$human
    humanGenes <- intersect(rownames(ref), human_mouse_gene_symbols$human)
    ref <- ref[humanGenes, ]
    rownames(ref) <- human_mouse_gene_symbols[humanGenes, ]$mouse
  }
  
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
MakeQuantiles <- function(ref, labels, probs, numThreads) {
  param <- BiocParallel::MulticoreParam(workers = numThreads)
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
  }, BPPARAM = param)
  names(quantilesMatList) <- cellTypes
  
  return(quantilesMatList)
}

# Function to create signatures for each cell type
CreateSignatures <- function(labels, depList, quantilesMatrix, probs, diffVals, minGenes,
                             maxGenes, minFracCtPassed, numThreads) {
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
  
  param <- BiocParallel::MulticoreParam(workers = numThreads)
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
  }, BPPARAM = param)
  
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
LearnParams <- function(gepMat, corMat, signatures, depList, topSpillValue, numThreads) {
  param <- BiocParallel::MulticoreParam(workers = numThreads)
  
  cellTypes <- colnames(gepMat)
  gepMatLinear <- 2^gepMat
  simFracs <- c(0, seq(0.01, 0.25, 0.01))
  fracToUse <- 0.25
  
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
  }, BPPARAM = param)
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
  }, BPPARAM = param)
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
  }, BPPARAM = param)
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


#' xCell2Train function
#'
#' This function generates a custom xCell2 reference object for cell type enrichment analysis.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import readr
#' @import BiocParallel
#' @importFrom utils data
#' @importFrom SummarizedExperiment assays colData
#' @importFrom minpack.lm nlsLM
#' @importFrom Rfast rowMedians rowmeans rowsums Sort
#' @importFrom Matrix rowMeans rowSums colSums Diagonal
#' @importFrom singscore rankGenes simpleScore
#' @importFrom stats coef cor lm quantile var
#' @importFrom methods new
#' @param ref A reference gene expression matrix (with genes in rows and samples/cells in columns),
#'        a SummarizedExperiment object, or a SingleCellExperiment object containing the expression
#'        data and sample metadata. If a SummarizedExperiment or SingleCellExperiment object is provided,
#'        the expression matrix should be stored in the "counts" slot of the `assays` component, 
#'        and the sample metadata (equivalent to the "labels" parameter) should be stored in `colData`.
#' @param mix A bulk mixture of gene expression data (genes in rows, samples in columns) (optional).
#' @param labels A data frame in which the rows correspond to samples in the ref.
#'  The data frame must have four columns:
#'  "ont": the cell type ontology as a character (i.e., "CL:0000545" or NA if there is no ontology).
#'  "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'  "sample": the cell type sample/cell that match the column name in ref.
#'  "dataset": sample's source dataset or subject (can be the same for all samples if no such information).
#'  This parameter is not needed if ref is a SummarizedExperiment or SingleCellExperiment object,
#'  as it should already be included in `colData`.
#' @param refType The reference gene expression data type: "rnaseq" for bulk RNA-Seq, "array" for micro-array, or "sc" for scRNA-Seq.
#' @param minPbCells For scRNA-Seq reference only - minimum number of cells in the pseudo-bulk (optional, default: 30).
#' @param minPbSamples For scRNA-Seq reference only - minimum number of pseudo-bulk samples (optional, default: 10).
#' @param minScGenes For scRNA-Seq reference only - minimum number of genes for pseudo-bulk samples (default: 10000).
#' @param useOntology A Boolean for considering cell type dependencies by using ontological integration (default: TRUE).
#' @param lineageFile Path to the cell type lineage file generated with `xCell2GetLineage` function and reviewed manually (optional).
#' @param numThreads Number of threads for parallel processing (default: 1).
#' @param humanToMouse A Boolean for converting human genes to mouse genes (default: FALSE).
#' @param topSpillValue Maximum spillover compensation correction value (default: 0.5).
#' @param returnSignatures A Boolean to return just the signatures (default: FALSE).
#' @param returnAnalysis A Boolean to return the xCell2Analysis results (do not return reference object) (default: FALSE).
#' @param useSpillover A Boolean to use spillover correction in xCell2Analysis (returnAnalysis must be TRUE) (default: TRUE).
#' @param spilloverAlpha A numeric for spillover alpha value in xCell2Analysis (returnAnalysis must be TRUE) (default: 0.5).
#' @return An S4 object containing cell types' signatures, linear transformation parameters, spillover matrix and dependencies.
#' @examples
#' # For detailed example read xCell2 vignette.
#'
#' # Extract reference matrix
#' data(dice_demo_ref, package = "xCell2")
#' dice_ref <- as.matrix(dice_demo_ref@assays@data$logcounts)
#' colnames(dice_ref) <- make.unique(colnames(dice_ref)) # Make samples samples unique
#'
#' # Extract reference metadata
#' dice_labels <- as.data.frame(dice_demo_ref@colData)
#'
#' # Prepare labels data frame
#' dice_labels$ont <- NA
#' dice_labels$sample <- colnames(dice_ref)
#' dice_labels$dataset <- "DICE"
#'
#' # Assign cell type ontology (optional but recommended)
#' dice_labels[dice_labels$label == "B cells", ]$ont <- "CL:0000236"
#' dice_labels[dice_labels$label == "Monocytes", ]$ont <- "CL:0000576"
#' dice_labels[dice_labels$label == "NK cells", ]$ont <- "CL:0000623"
#' dice_labels[dice_labels$label == "T cells, CD8+", ]$ont <- "CL:0000625"
#' dice_labels[dice_labels$label == "T cells, CD4+", ]$ont <- "CL:0000624"
#' dice_labels[dice_labels$label == "T cells, CD4+, memory", ]$ont <- "CL:0000897"
#'
#' # Reproducibility (optional): Set seed before running `xCell2Train`  as generating pseudo-bulk
#' # samples from scRNA-Seq reference based on random sampling of cells.
#' set.seed(123)
#'
#' # Generate custom xCell2 reference object
#' DICE.xCell2Ref <- xCell2::xCell2Train(ref = dice_ref, labels = dice_labels, refType = "rnaseq")
#'
#' @export
xCell2Train <- function(ref,
                        mix = NULL,
                        labels = NULL,
                        refType,
                        humanToMouse = FALSE,
                        lineageFile = NULL,
                        numThreads = 1,
                        useOntology = TRUE,
                        returnSignatures = FALSE,
                        returnAnalysis = FALSE,
                        useSpillover = TRUE,
                        spilloverAlpha = 0.5,
                        minPbCells = 30,
                        minPbSamples = 10,
                        minScGenes = 1e4,
                        topSpillValue = 0.5) {
  
  
  inputsValidated <- ValidateInputs(ref, labels, refType)
  ref <- inputsValidated$ref
  labels <- inputsValidated$labels
  
  # Generate pseudo bulk from scRNA-Seq reference
  if (refType == "sc") {
    message("Converting scRNA-seq reference to pseudo bulk...")
    pbData <- ScToPseudoBulk(ref, labels, minPbCells, minPbSamples)
    ref <- pbData$ref
    labels <- pbData$labels
  }
  
  # Prepare reference and mixture data
  out <- PrepRefMix(ref, mix, refType, minScGenes, humanToMouse)
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
  quantilesMatrix <- MakeQuantiles(ref, labels, probs, numThreads)
  signatures <- CreateSignatures(
    labels, depList, quantilesMatrix, probs, diffVals, minGenes,
    maxGenes, minFracCtPassed, numThreads
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
  params <- LearnParams(gepMat, corMat, signatures, depList, topSpillValue, numThreads)
  
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
                                  spilloverAlpha = spilloverAlpha, numThreads = numThreads
    )
    return(res)
  } else {
    return(xCell2S4)
  }
}
