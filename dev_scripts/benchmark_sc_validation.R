library(babelwhale)
library(rhdf5)
library(parallel)
library(RcppHungarian)
library(stringr)

setwd("/bigdata/almogangel/twelve_years_decon_paper/analysis/")
source('./include/utils.R')


PATH_TO_MATLAB_LICENSE <- "/path/to/matlab/license.lic"
DOCKER_TAG_PREFIX <- "deconvolution/"
DOCKER_TAG_SUFFIX <- ":latest"

methods <- c("xCell2", "EPIC", "MuSic", "DeconRNASeq", "dtangle", "scaden", "quanTIseq", "MCPcounter")

dataFiles <- list.files('./data/sim')
allConfig <- expand.grid(dataFiles, methods)

# 8 cores, 16GB per job, 12GB softlimit, 12 concurrent jobs, max runtime of 12 hours

dir.create('./results/accuracy', showWarnings = FALSE, recursive = TRUE)
lapply(seq_len(nrow(allConfig)), function(i) {
  dataFile <- paste0('./data/sim/', as.character(allConfig[i,]$Var1))
  method <- as.character(allConfig[i,]$Var2)
  # method <- "xCell2"

  resFile <- paste0('./results/accuracy/', method, '_', as.character(allConfig[i,]$Var1))
  if (file.exists(resFile)) {
    return(NULL)
  }

  data <- readRDS(dataFile)
  groundTruth <- data$bulkRatio
  data <- data[names(data) != "bulkRatio"]

  start <- Sys.time()
  res <- do.call(
    runDeconvolution,
    c(
      list(methods = method, verbose = TRUE,
           dockerArgs = c(
             '--cpus=8.0',
             '-m=16G',
             '--memory-reservation=12G'
           ),
           timeout = 12*3600,
           matlabLicenseFile=PATH_TO_MATLAB_LICENSE),
      data
    )
  )
  runningTime <- Sys.time() - start

  res <- res[[method]]
  res$groundTruth <- groundTruth
  res$runningTime <- runningTime

  saveRDS(res, file = resFile)

  NULL
})

# gather results

dataPath <- "./data/sim"
datasets <- list.files(dataPath) %>%
  str_match('([^/]+).rds$') %>%
  { function(x) x[, 2] }()

methods <- ALL_METHODS
rawRes <- lapply(datasets, function(dataset) {
  groundTruthFile <- paste0('./data/sim/', dataset, '.rds')
  groundTruth     <- readRDS(groundTruthFile)$bulkRatio %>% t()
  colnames(groundTruth) <- gsub("_", "-", colnames(groundTruth)) # Almog


  res <- lapply(methods, function(method) {
    # message(dataset, method)
    resFile <- paste0('./results/accuracy/', method, '_', dataset, '.rds')
    if (!file.exists(resFile)) return(NULL)
    print(resFile)

    res <- readRDS(resFile)

    if (!is.null(res$stderr)) return(NULL)

    if (!rownames(res$groundTruth)[1] %in% colnames(res$P)) {
      if (ncol(res$P) < nrow(res$groundTruth))
      {
        return(NULL)
      }
      if (is.null(colnames(res$P))) {
        colnames(res$P) <- paste0('CT', seq_len(ncol(res$P)))
      }
    }

    P <- res$P
    P <- apply(P, 1, function(x) {
      x <- x - min(x)
      x <- x / sum(x)
      x
    }) %>% t()

    colnames(P) <- gsub("_", "-", colnames(P)) # Almog

    if (!colnames(groundTruth)[1] %in% colnames(P)) {
      mapRes      <- detectCellType(P, t(groundTruth))
      colnames(P) <- mapRes[colnames(P),]$to
    }

    list(
      P           = P[, colnames(groundTruth)] %>%
        data.frame(
          method  = method,
          dataset = dataset,
          sample  = rownames(P)
        ) %>%
        tidyr::gather("cellType", "value", -method, -dataset, -sample),
      groundTruth = groundTruth %>%
        data.frame(
          method  = method,
          dataset = dataset,
          sample  = rownames(groundTruth)
        ) %>%
        tidyr::gather("cellType", "value", -method, -dataset, -sample)
    )
  })

  list(
    P           = lapply(res, function(x) x$P) %>% do.call(what = rbind),
    groundTruth = lapply(res, function(x) x$groundTruth) %>% do.call(what = rbind)
  )
}) %>% { function(res) {
  list(
    P           = lapply(res, function(x) x$P) %>% do.call(what = rbind),
    groundTruth = lapply(res, function(x) x$groundTruth) %>% do.call(what = rbind)
  )
} }()

datasets <- datasets[datasets != "Liver"] # TODO: fix liver

allRes <- lapply(datasets, function(dataset) {
  print(dataset)
  groundTruthFile <- paste0('./data/sim/', dataset, '.rds')
  groundTruth     <- readRDS(groundTruthFile)$bulkRatio %>% t()



  # method <- methods[1]
  lapply(methods, function(method) {
    resFile <- paste0('./results/accuracy/', method, '_', dataset, '.rds')
    if (!file.exists(resFile)) return(NULL)

    res <- readRDS(resFile)

    colnames(res$P) <- gsub("-", "_", colnames(res$P)) # Almog


    if (!is.null(res$stderr)) return(NULL)

    if (!rownames(res$groundTruth)[1] %in% colnames(res$P)) {
      if (ncol(res$P) < nrow(res$groundTruth))
      {
        return(NULL)
      }
      if (is.null(colnames(res$P))) {
        colnames(res$P) <- paste0('CT', seq_len(ncol(res$P)))
      }
    }

    SCorr <- mean(evaluateResult(res$P, res$groundTruth, method = "sample"), na.rm = TRUE)
    CCorr <- mean(evaluateResult(res$P, res$groundTruth, method = "CT"), na.rm = TRUE)
    MAE   <- mean(evaluateResult(res$P, res$groundTruth, metric = "MAE"), na.rm = TRUE)
    MAECorr   <- evaluateResult(res$P, res$groundTruth, metric = "samplePairwise")

    data.frame(
      dataset = dataset,
      method  = method,
      SCorr   = SCorr,
      CCorr   = CCorr,
      MAE     = MAE,
      MAECorr     = MAECorr
    )
  }) %>% do.call(what = rbind)

}) %>% do.call(what = rbind)

saveRDS(allRes, file = './results/accuracy/allRes.rds')
