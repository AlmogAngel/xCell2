library(DeconBenchmark)
library(tidyverse)

makeSigMatForCIBERSORTx <- function(singleCellExpr, refName){
  dir <- "/bigdata/almogangel/CIBERSORTx_docker"

  # Make pheno file
  all_celltypes <- unique(colnames(singleCellExpr))
  pheno.mat <- t(sapply(all_celltypes, function(ct){
    ifelse(ct == colnames(singleCellExpr), 1, 2)
  }))
  pheno.df <- data.frame(cbind(rownames(pheno.mat), pheno.mat))
  pheno_file <- paste0(dir, "/", refName, "_pheno.txt")

  if(!file.exists(pheno_file)){
    write.table(pheno.df, file = pheno_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }

  # Make reference file
  tmp <- cbind("genes" = rownames(singleCellExpr), singleCellExpr)
  tmp <- singleCellExpr %>%
    as.data.frame(.) %>%
    rownames_to_column("genes")
  refsample_file <- paste0(dir, "/", refName, "_ref.txt")
  if(!file.exists(refsample_file)){
    write.table(tmp, file = refsample_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }

  # Make signature matrix
  token <- "b72da36961922443b75a1b65beef27c0"

  # Make results directory
  results_dir <- paste0(dir, "/results")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm -f ", results_dir, "/*"))
  }

  cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                token, " --single_cell FALSE --refsample ", refsample_file, " --phenoclasses ", pheno_file, " --QN TRUE 1> ",
                results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)

  sigmat_path <- paste0(dir, "/", refName, "_sigmat.txt")
  file.copy(paste0(results_dir, "/CIBERSORTx_", refName, "_pheno.CIBERSORTx_", refName, "_ref.bm.K999.txt"), sigmat_path)

  return(sigmat_path)
}

runCIBERSORTx <- function(mix, sigmat_path, dir = "/bigdata/almogangel/CIBERSORTx_docker"){

  token <- "b72da36961922443b75a1b65beef27c0"

  mix_tmp <- cbind("genes" = rownames(mix), mix)
  mix_file <- paste0(dir, "/mix-tmp.txt")
  write.table(mix_tmp, file = mix_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

  # Make results directory
  results_dir <- paste0(dir, "/results")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm -f ", results_dir, "/*"))
  }

  cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                token, " --single_cell FALSE --sigmatrix ", sigmat_path, " --mixture ", mix_file, " --QN TRUE 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)

  # Load results
  cibersortx_out <- t(read.table(paste0(results_dir,  "/CIBERSORTx_Results.txt"), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
  cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]


  return(cibersortx_out)
}

runMCPcounter <- function(expression, markers){

  appendSignatures=function(xp,markers){
    res=as.data.frame(do.call(cbind,
                              lapply(markers,function(x){
                                apply(xp[intersect(row.names(xp),x),,drop=F],2,mean,na.rm=T)
                              })))
    res
  }


  t(appendSignatures(expression,markers))

}

generateReferenceMarkersBulk  <- function(singleCellExpr, singleCellLabels, log2Threshold = 1) {

  singleCellLabels <- as.character(singleCellLabels)
  colnames(singleCellExpr) <- singleCellLabels

  #for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
  keeps <- sapply(unique(singleCellLabels), function(cellType) {
    hits <- singleCellLabels == cellType
    rowSums(singleCellExpr[, hits, drop = FALSE] != 0) >= ceiling(0.3 * sum(hits))
  })

  #normalization
  singleCellExpr <- singleCellExpr[rowSums(keeps) > 0,]
  singleCellExpr <- edgeR::DGEList(singleCellExpr)
  singleCellExpr <- edgeR::calcNormFactors(singleCellExpr, method = "TMM")

  # INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account
  #[compare one group with average expression of all other groups]
  annotation <- factor(singleCellLabels)
  design <- model.matrix(~0 + annotation)
  colnames(design) <- sapply(strsplit(colnames(design), "annotation"), function(x) x[2])
  constrastMatrix <- matrix((-1 / ncol(design)), nrow = ncol(design), ncol = ncol(design))
  colnames(constrastMatrix) <- colnames(design)
  diag(constrastMatrix) <- (ncol(design) - 1) / ncol(design)

  v <- limma::voom(singleCellExpr, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::contrasts.fit(fit, constrastMatrix)
  fit <- limma::eBayes(fit, trend = TRUE)


  topTableResults <- limma::topTable(fit, coef = seq_len(ncol(constrastMatrix)), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2Threshold)
  topTableResults <- topTableResults[, 1:(ncol(topTableResults) - 4)]

  ERCCGenes <- grep("ERCC-", rownames(topTableResults))
  if (length(ERCCGenes) > 0) {
    topTableResults <- topTableResults[-ERCCGenes,]
  }

  markers <- apply(topTableResults, 1, function(x) {
    temp <- sort(x)
    ((temp[ncol(topTableResults)] - temp[ncol(topTableResults) - 1]) >= log2Threshold) |
      (abs(temp[1] - temp[2]) >= log2Threshold)
  })

  topTableResults <- topTableResults[markers,]

  markers <- cbind.data.frame(
    rownames(topTableResults),
    t(apply(topTableResults, 1, function(x) {
      temp <- max(x)
      if (temp < log2Threshold) {
        temp <- c(min(x), colnames(topTableResults)[which.min(x)])
      } else {
        temp <- c(max(x), colnames(topTableResults)[which.max(x)])
      }
      temp
    }))
  )

  colnames(markers) <- c("gene", "log2FC", "cellType")
  genes <- as.character(markers$gene)
  cellTypes <- as.character(markers$cellType)

  markers <- lapply(unique(cellTypes), function(cellType) genes[cellTypes == cellType])
  names(markers) <- unique(cellTypes)

  markers
}

generateReferenceBulk <- function(singleCellExpr, singleCellLabels, types = c("markers", "sigGenes", "signature", "cellTypeExpr"), log2Threshold = 1) {
  reference <- list()

  if (any(c("markers", "sigGenes", "signature") %in% types)) {
    reference$markers <- generateReferenceMarkersBulk(singleCellExpr, singleCellLabels, log2Threshold)
  }

  if (any(c("sigGenes", "signature") %in% types)) {
    reference$sigGenes <- unique(unlist(reference$markers))
  }

  if (any(c("signature", "cellTypeExpr") %in% types)) {
    cellTypeExpr <- matrix(ncol = length(unique(singleCellLabels)), nrow = nrow(singleCellExpr))
    colnames(cellTypeExpr) <- unique(singleCellLabels)
    rownames(cellTypeExpr) <- rownames(singleCellExpr)

    for (cellType in unique(singleCellLabels)) {
      #tmp <- rowSums(singleCellExpr[, singleCellLabels == cellType])
      #tmp <- tmp / sum(tmp) * 1e6
      tmp <- apply(singleCellExpr[, singleCellLabels == cellType], 1, median)
      cellTypeExpr[, cellType] <- tmp
    }

    reference$cellTypeExpr <- cellTypeExpr
  }

  if ("signature" %in% types) {
    reference$signature <- reference$cellTypeExpr[reference$sigGenes, ]
  }

  reference[types]
}


# Load reference
labels <- readRDS("/bigdata/almogangel/xCell2/dev_data/sref_blood_labels_bulk.rds")
ref <- readRDS("/bigdata/almogangel/xCell2/dev_data/sref_blood_data_bulk.rds")

# Subset cell types (must be independent cell types for deconvolution methods)
celltypes2use <- c("B cell", "CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell",
                   "basophil", "monocyte", "natural killer cell", "neutrophil", "eosinophil")
sample2keep <- labels$label %in% celltypes2use
labels <- labels[sample2keep,]
ref <- ref[,sample2keep]

# Change labels
labels$label <- gsub("-", "_", labels$label)
labels$label <- gsub(" ", "_", labels$label)
labels$label <- gsub(",", "_", labels$label)

singleCellExpr <- ref
singleCellLabels <- labels$label
colnames(singleCellExpr) <- singleCellLabels
singleCellSubjects <- labels$dataset

# Make signature matrix for CIBERSORTx
# sigmat_path <- makeSigMatForCIBERSORTx(singleCellExpr, refName = "sref_blood")
sigmat_path <- "/bigdata/almogangel/CIBERSORTx_docker/sref_blood_sigmat.txt"

# Make reference for other methods
# reference <- generateReferenceBulk(singleCellExpr, singleCellLabels, types = c("markers", "sigGenes", "signature", "cellTypeExpr"))
# saveRDS(reference, "/bigdata/almogangel/xCell2/dev_data/DeconBenchmark_sref_blood_reference.rds")
reference <- readRDS("/bigdata/almogangel/xCell2/dev_data/DeconBenchmark_sref_blood_reference.rds")

truths_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/"
mix_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/"
blood_ds <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813")

results.list <- list()
for (file in blood_ds) {

  print(file)

  bulk <- as.matrix(read.table(paste0(mix_dir, file, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))

  # Run EPIC
  epic.ref <- list(
    refProfiles = reference$cellTypeExpr,
    sigGenes = reference$sigGenes
    # TODO: add refProfiles.var
  )
  epic.out <- EPIC::EPIC(bulk, epic.ref)$mRNAProportions
  results.list[[file]][["EPIC"]] <- epic.out

  # Run DeconRNASeq
  DeconRNASeq.out <- runDeconvolution("DeconRNASeq", bulk = bulk, signature = reference$signature)
  results.list[[file]][["DeconRNASeq"]] <- DeconRNASeq.out$DeconRNASeq$P

  # Run dtangle
  genes2use <- intersect(rownames(bulk), rownames(reference$cellTypeExpr))
  dtangle.out <- dtangle::dtangle(Y=t(bulk[genes2use,]), references = t(reference$cellTypeExpr[genes2use,]))$estimates
  results.list[[file]][["dtangle"]] <- dtangle.out

  # Run scaden
  # TODO: runtime too long
  scaden.out <- runDeconvolution("scaden", bulk = bulk, singleCellExpr = singleCellExpr, singleCellLabels = singleCellLabels)
  results.list[[file]][["scaden"]] <- scaden.out

  # Run MCPcounter
  MCPcounter.out <- t(runMCPcounter(bulk, reference$markers))
  results.list[[file]][["MCPcounter"]] <- MCPcounter.out

  # Run CIBERSORTx
  CIBERSORTx.out <- runCIBERSORTx(mix = bulk, sigmat_path = sigmat_path, dir = "/bigdata/almogangel/CIBERSORTx_docker")
  results.list[[file]][["CIBERSORTx"]] <- t(CIBERSORTx.out)


}

saveRDS(results.list, paste0("/bigdata/almogangel/twelve_years_decon_paper/sref_blood_deconvolution_results.rds"))



