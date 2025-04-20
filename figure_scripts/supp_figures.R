library(tidyverse)
library(GEOquery)
library(R.utils)
library(Biobase)
library(xCell2)
library(BiocParallel)

options(timeout = 1200) 

source("/bigdata/almogangel/xCell2_dev/paper_figures/load_benchmark_data.R")

# Robustness to normalization -------------------

geneids2symbols <- function(expr_data){
  
  if (any(duplicated(rownames(expr_data)))) {
    errorCondition("Duplicated gene IDs!")
    return(NA)
  }
  
  shared_ids <- intersect(rownames(expr_data), rownames(gene_ann))
  expr_data <- expr_data[shared_ids,]
  rownames(expr_data) <- gene_ann[shared_ids,]$Symbol
  
  return(expr_data)
}
download_and_load_gse <- function(url = NULL, gse_id, dest_dir = "/bigdata/almogangel/xCell2_data/robust_norm_analysis/", annotate_genes = FALSE) {
  
  if (is.null(url)) {
    dir <- paste0(getwd(), "/data/", gse_id)
    if (!dir.exists(dir)) {
      dir.create(dir)
      gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = dir)
    }else{
      gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = dir)
    }
    return(gset)
  }
  
  
  # Create a directory named after the GSE ID
  output_dir <- file.path(dest_dir, gse_id)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Define the output file path
  output_file <- file.path(output_dir, paste0(gse_id, ".gz"))
  unzipped_file <- sub("\\.gz$", "", output_file)  # Remove the .gz extension
  unzipped_file <- paste0(unzipped_file, ".tsv")
  
  # Download the file
  if (!file.exists(unzipped_file)) {
    download.file(url, destfile = output_file, mode = "wb")
    message("File downloaded: ", output_file)
    
    # Gunzip the file
    gunzip(output_file, destname = unzipped_file, overwrite = TRUE)
    message("File unzipped: ", unzipped_file)
  } 
  
  # Read the gene expression matrix into R
  expr_data <- read.delim(unzipped_file, header = TRUE, sep = "\t", row.names = 1)
  if (ncol(expr_data) == 0) {
    expr_data <- read.delim(unzipped_file, header = TRUE, sep = ",", row.names = 1)
  }
  
  if (annotate_genes) {
    if (all(startsWith(prefix = "ENSG", rownames(expr_data)))) {
      expr_data <- convert_ensembl_to_symbol(expr_data)
    }else{
      expr_data <- geneids2symbols(expr_data)
    }
  }
  
  # Check the data and return
  message("Expression matrix loaded. Dimensions: ", dim(expr_data)[1], " x ", dim(expr_data)[2])
  return(expr_data)
}
get_expression_matrix <- function(accession, dest_dir = "/bigdata/almogangel/xCell2_data/robust_norm_analysis/", normalization = "RMA") {
  # Check required packages
  required_pkgs <- c("GEOquery", "affy", "Biobase")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Please install the '", pkg, "' package.")
    }
  }
  
  # Create a directory for the accession if it doesn't exist
  acc_dir <- file.path(dest_dir, accession)
  if (!dir.exists(acc_dir)) {
    dir.create(acc_dir, recursive = TRUE)
  }
  
  # Check if raw files already exist (either a *_RAW.tar or any CEL files)
  raw_tar_files <- list.files(acc_dir, pattern = "_RAW\\.tar$", full.names = TRUE)
  cel_files   <- list.files(acc_dir, pattern = "\\.CEL$|\\.cel$|\\.CEL\\.gz$|\\.cel\\.gz$", full.names = TRUE)
  
  if (length(raw_tar_files) == 0 && length(cel_files) == 0) {
    message("No raw files found in ", acc_dir, ". Downloading supplemental files for ", accession, "...")
    GEOquery::getGEOSuppFiles(accession, baseDir = dest_dir)
    raw_tar_files <- list.files(acc_dir, pattern = "_RAW\\.tar$", full.names = TRUE)
    if (length(raw_tar_files) == 0) {
      stop("No *_RAW.tar file found for accession ", accession)
    }
  } else {
    message("Raw files already exist in ", acc_dir, ". Skipping download.")
  }
  
  # If a *_RAW.tar file exists, extract it only if no CEL files are found yet
  if (length(raw_tar_files) > 0 && length(cel_files) == 0) {
    tar_file <- raw_tar_files[1]
    message("Extracting CEL files from ", basename(tar_file), " ...")
    utils::untar(tar_file, exdir = acc_dir)
  } else {
    message("CEL files already extracted (or available).")
  }
  
  
  # Find CEL files including gzipped ones
  cel_files <- list.files(acc_dir, pattern = "\\.CEL$|\\.cel$|\\.CEL\\.gz$|\\.cel\\.gz$", full.names = TRUE)
  if (length(cel_files) == 0) {
    stop("No CEL files found in the extracted tar archive.")
  }
  
  # Decompress any gzipped CEL files
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Please install the 'R.utils' package.")
  }
  for (f in cel_files) {
    if (grepl("\\.gz$", f)) {
      message("Decompressing file: ", basename(f))
      out_file <- sub("\\.gz$", "", f)
      if (!file.exists(out_file)) {
        R.utils::gunzip(f, destname = out_file, overwrite = TRUE)
      }
    }
  }
  
  # Update list of CEL files to include only the decompressed ones
  cel_files <- list.files(acc_dir, pattern = "\\.CEL$|\\.cel$", full.names = TRUE)
  
  message("Reading CEL files...")
  affy_data <- affy::ReadAffy(filenames = cel_files)
  
  # Normalize using the selected method
  if (normalization == "RMA") {
    message("Normalizing data using RMA...")
    eset <- affy::rma(affy_data)
  } else if (normalization == "MAS5") {
    message("Normalizing data using MAS5...")
    eset <- affy::mas5(affy_data)
  } else if (normalization == "QN") {
    message("Normalizing data using only Quantile Normalization...")
    # Use expresso() with background correction and PM correction turned off
    eset <- affy::expresso(affy_data, 
                           bg.correct = FALSE,
                           normalize = TRUE,
                           normalize.method = "quantiles",
                           pmcorrect.method = "pmonly",
                           summary.method = "medianpolish")
  } else {
    stop("Normalization method not supported. Choose RMA, MAS5, or QN.")
  }
  
  # Extract the normalized expression matrix
  expr_matrix <- Biobase::exprs(eset)
  
  # Convert probe annotations to gene symbols if possible.
  # This example uses the hgu133plus2 platform. If the platform is different,
  # adjust the annotation package accordingly.
  platform <- Biobase::annotation(eset)
  if (!is.null(platform) && platform == "hgu133plus2") {
    if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
      warning("hgu133plus2.db package not installed; returning data with probe IDs.")
    } else {
      message("Mapping probe IDs to gene symbols using hgu133plus2.db...")
      # Map probe IDs to gene symbols
      probe_ids <- rownames(expr_matrix)
      gene_symbols <- AnnotationDbi::mapIds(hgu133plus2.db::hgu133plus2.db,
                                            keys = probe_ids,
                                            column = "SYMBOL",
                                            keytype = "PROBEID",
                                            multiVals = "first")
      # Filter out probes with no mapping
      valid <- !is.na(gene_symbols)
      expr_matrix <- expr_matrix[valid, , drop = FALSE]
      gene_symbols <- gene_symbols[valid]
      
      # Create a data frame, add gene symbol column, and aggregate by gene symbol
      df <- as.data.frame(expr_matrix)
      df$gene <- gene_symbols
      agg_df <- aggregate(. ~ gene, data = df, FUN = mean)
      rownames(agg_df) <- agg_df$gene
      agg_df$gene <- NULL
      expr_matrix <- as.matrix(agg_df)
    }
  } else {
    message("No recognized annotation available for mapping probes to gene symbols; returning data with original probe IDs.")
  }
  

  return(expr_matrix)
}

gene_ann <- read.delim("/bigdata/almogangel/tme_treatment_response/Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t")
gene_ann <- gene_ann[gene_ann$GeneType == "protein-coding",]
rownames(gene_ann) <- as.character(gene_ann$GeneID)

human.refs <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/human_references.rds")
human.vals <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/human_validation.rds")
refs.vals.matched.human <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/refs.vals.matched.human.rds")

human.vals.tpm <- human.vals
human.vals.counts <- human.vals
human.vals.fpkm <- human.vals

# RNA-Seq

# GSE107011
GSE107011.tpm  <- human.vals.tpm$mixtures$blood$GSE107011
colSums(GSE107011.tpm)
dim(GSE107011.tpm)

GSE107011.counts <- download_and_load_gse(url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE107011&format=file&file=GSE107011_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                                          gse_id = "GSE107011-counts", annotate_genes = TRUE)
GSE107011.fpkm <- download_and_load_gse(url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE107011&format=file&file=GSE107011_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz",
                                        gse_id = "GSE107011-fpkm", annotate_genes = TRUE)
colSums(GSE107011.counts)
colSums(GSE107011.fpkm)

gset <- getGEO("GSE107011", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "/bigdata/almogangel/xCell2_data/robust_norm_analysis/GSE107011-counts/")
m <- pData(gset$GSE107011_series_matrix.txt.gz)
m <- m[grepl(pattern = "PBMC", m$title),]
rownames(m) <- gsub("_.*", "", m$title)
m <- m[colnames(human.vals.tpm$mixtures$blood$GSE106898),]


GSE107011.counts <- GSE107011.counts[,m$geo_accession]
colnames(GSE107011.counts) <- colnames(GSE107011.tpm)
colSums(GSE107011.counts)
dim(GSE107011.counts)

GSE107011.fpkm <- GSE107011.fpkm[,m$geo_accession]
colnames(GSE107011.fpkm) <- colnames(GSE107011.tpm)
colSums(GSE107011.fpkm)
dim(GSE107011.fpkm)

genes2use <- intersect(rownames(GSE107011.tpm), rownames(GSE107011.counts))
GSE107011.tpm <- GSE107011.tpm[genes2use,]
GSE107011.counts <- GSE107011.counts[genes2use,]
GSE107011.counts <- as.matrix(GSE107011.counts)
GSE107011.fpkm <- GSE107011.fpkm[genes2use,]
GSE107011.fpkm <- as.matrix(GSE107011.fpkm)
dim(GSE107011.tpm)
dim(GSE107011.counts)
dim(GSE107011.fpkm)

human.vals.tpm$mixtures$blood$GSE107011 <- GSE107011.tpm
human.vals.counts$mixtures$blood$GSE107011 <- GSE107011.counts
human.vals.fpkm$mixtures$blood$GSE107011 <- GSE107011.fpkm


# GSE130824
GSE130824.tpm  <- human.vals.tpm$mixtures$blood$GSE130824
colSums(GSE130824.tpm)
dim(GSE130824.tpm)
GSE130824.counts <- download_and_load_gse(url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE130824&format=file&file=GSE130824_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                                          gse_id = "GSE130824-counts", annotate_genes = TRUE)
colSums(GSE130824.counts)
dim(GSE130824.counts)
GSE130824.fpkm <- download_and_load_gse(url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE130824&format=file&file=GSE130824_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz",
                                        gse_id = "GSE130824-fpkm", annotate_genes = TRUE)
colSums(GSE130824.fpkm)
dim(GSE130824.fpkm)

gset <- getGEO("GSE130824", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "/bigdata/almogangel/xCell2_data/robust_norm_analysis/GSE130824-counts/")
m <- pData(gset$GSE130824_series_matrix.txt.gz)
samples <- m$title
samples <- gsub(":.*", "", samples)
samples <- gsub("-", "", samples)
samples <- gsub("B0", "BD", samples)

rownames(m) <- samples
m <- m[colnames(GSE130824.tpm),]

GSE130824.counts <- GSE130824.counts[,m$geo_accession]
colnames(GSE130824.counts) <- colnames(GSE130824.tpm)
colSums(GSE130824.counts)
dim(GSE130824.counts)

GSE130824.fpkm <- GSE130824.fpkm[,m$geo_accession]
colnames(GSE130824.fpkm) <- colnames(GSE130824.tpm)
colSums(GSE130824.fpkm)
dim(GSE130824.fpkm)

colSums(GSE130824.tpm)
dim(GSE130824.tpm)

genes2use <- intersect(rownames(GSE130824.tpm), rownames(GSE130824.counts))
GSE130824.tpm <- GSE130824.tpm[genes2use,]
GSE130824.counts <- GSE130824.counts[genes2use,]
GSE130824.counts <- as.matrix(GSE130824.counts)
GSE130824.fpkm <- GSE130824.fpkm[genes2use,]
GSE130824.fpkm <- as.matrix(GSE130824.fpkm)
dim(GSE130824.tpm)
dim(GSE130824.counts)
dim(GSE130824.fpkm)

human.vals.tpm$mixtures$blood$GSE130824 <- GSE130824.tpm
human.vals.counts$mixtures$blood$GSE130824 <- GSE130824.counts
human.vals.fpkm$mixtures$blood$GSE130824 <- GSE130824.fpkm

# GSE64655
GSE64655.tpm  <- human.vals.tpm$mixtures$blood$GSE64655
colSums(GSE64655.tpm)
dim(GSE64655.tpm)
GSE64655.counts <- download_and_load_gse(url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE64655&format=file&file=GSE64655_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                                          gse_id = "GSE64655-counts", annotate_genes = TRUE)
colSums(GSE64655.counts)
dim(GSE64655.counts)
GSE64655.fpkm <- download_and_load_gse(url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE64655&format=file&file=GSE64655_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz",
                                        gse_id = "GSE64655-fpkm", annotate_genes = TRUE)
colSums(GSE64655.fpkm)
dim(GSE64655.fpkm)

gset <- getGEO("GSE64655", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "/bigdata/almogangel/xCell2_data/robust_norm_analysis/GSE64655-counts/")
m <- pData(gset$GSE64655_series_matrix.txt.gz)
samples <- gsub("_SL.*", "", m$title)
samples <- gsub("PBMC_", "PBMC_d", samples)
samples <- gsub("HD30_PMBC_0", "HD30_PBMC_d0", samples)
rownames(m) <- samples
m <- m[colnames(GSE64655.tpm),]

GSE64655.counts <- GSE64655.counts[,m$geo_accession]
colnames(GSE64655.counts) <- colnames(GSE64655.tpm)
colSums(GSE64655.counts)
dim(GSE64655.counts)

GSE64655.fpkm <- GSE64655.fpkm[,m$geo_accession]
colnames(GSE64655.fpkm) <- colnames(GSE64655.tpm)
colSums(GSE64655.fpkm)
dim(GSE64655.fpkm)

colSums(GSE64655.tpm)
dim(GSE64655.tpm)

genes2use <- intersect(rownames(GSE64655.tpm), rownames(GSE64655.counts))
GSE64655.tpm <- GSE64655.tpm[genes2use,]
GSE64655.counts <- GSE64655.counts[genes2use,]
GSE64655.counts <- as.matrix(GSE64655.counts)
GSE64655.fpkm <- GSE64655.fpkm[genes2use,]
GSE64655.fpkm <- as.matrix(GSE64655.fpkm)
dim(GSE64655.tpm)
dim(GSE64655.counts)
dim(GSE64655.fpkm)

human.vals.tpm$mixtures$blood$GSE64655 <- GSE64655.tpm
human.vals.counts$mixtures$blood$GSE64655 <- GSE64655.counts
human.vals.fpkm$mixtures$blood$GSE64655 <- GSE64655.fpkm

# GSE107572
GSE107572.tpm  <- human.vals.tpm$mixtures$blood$GSE107572
colSums(GSE107572.tpm)
dim(GSE107572.tpm)
GSE107572.counts <- download_and_load_gse(url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE107572&format=file&file=GSE107572_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                                         gse_id = "GSE107572-counts", annotate_genes = TRUE)
colSums(GSE107572.counts)
dim(GSE107572.counts)
GSE107572.fpkm <- download_and_load_gse(url = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE107572&format=file&file=GSE107572_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz",
                                       gse_id = "GSE107572-fpkm", annotate_genes = TRUE)
colSums(GSE107572.fpkm)
dim(GSE107572.fpkm)

gset <- getGEO("GSE107572", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "/bigdata/almogangel/xCell2_data/robust_norm_analysis/GSE107572-counts/")
m <- pData(gset$GSE107572_series_matrix.txt.gz)
rownames(m) <- colnames(GSE107572.tpm)

colnames(GSE107572.counts) <- colnames(GSE107572.tpm)
colSums(GSE107572.counts)
dim(GSE107572.counts)

colnames(GSE107572.fpkm) <- colnames(GSE107572.tpm)
colSums(GSE107572.fpkm)
dim(GSE107572.fpkm)

colSums(GSE107572.tpm)
dim(GSE107572.tpm)

genes2use <- intersect(rownames(GSE107572.tpm), rownames(GSE107572.counts))
GSE107572.tpm <- GSE107572.tpm[genes2use,]
GSE107572.counts <- GSE107572.counts[genes2use,]
GSE107572.counts <- as.matrix(GSE107572.counts)
GSE107572.fpkm <- GSE107572.fpkm[genes2use,]
GSE107572.fpkm <- as.matrix(GSE107572.fpkm)
dim(GSE107572.tpm)
dim(GSE107572.counts)
dim(GSE107572.fpkm)

human.vals.tpm$mixtures$blood$GSE107572 <- GSE107572.tpm
human.vals.counts$mixtures$blood$GSE107572 <- GSE107572.counts
human.vals.fpkm$mixtures$blood$GSE107572 <- GSE107572.fpkm

saveRDS(human.vals.tpm, "/bigdata/almogangel/xCell2_data/robust_norm_analysis/human.vals.tpm.rds")
saveRDS(human.vals.counts, "/bigdata/almogangel/xCell2_data/robust_norm_analysis/human.vals.counts.rds")
saveRDS(human.vals.fpkm, "/bigdata/almogangel/xCell2_data/robust_norm_analysis/human.vals.fpkm.rds")

# Run xCell2

ds2check <- c("GSE107011", "GSE130824", "GSE64655", "GSE107572")
params2use <- list(
  num_threads = 30,
  min_pb_cells = 20,
  min_pb_samples = 10,
  min_sc_genes = 1e4,
  use_ontology = TRUE,
  return_signatures = FALSE,
  return_analysis = FALSE,
  use_sillover = TRUE,
  spillover_alpha = 0.5
)

refs.vals.matched.human.xcell2 <- refs.vals.matched.human |> 
  mutate(n_val_samples = ncol(human.vals$truth[[val_type]][[val_dataset]])) |> 
  mutate(method = "xCell2", .before = everything()) |> 
  ungroup() |> 
  filter(val_dataset %in% ds2check)


dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/"
xcell2_tpm_res <- paste0(dir2use, "/fig1supp_xcell2_tpm_res.rds")
get_xcell2_benchmark_results(cyto.vals = human.vals.tpm, benchmark_table = refs.vals.matched.human.xcell2, save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 6, output_name = xcell2_tpm_res)
xcell2_tpm_res <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp//fig1supp_xcell2_tpm_res.rds")
benchmark_correlations.tpm <- get_xcell2_correlations(benchmark_table = refs.vals.matched.human.xcell2, xCell2results = xcell2_tpm_res, round_results = 3, just_xcell2 = TRUE, weight_cors = FALSE, by_val = FALSE)

dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/counts"
xcell2_counts_res <- paste0(dir2use, "/fig1supp_xcell2_counts_res.rds")
get_xcell2_benchmark_results(cyto.vals = human.vals.counts, benchmark_table = refs.vals.matched.human.xcell2, save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 6, output_name = xcell2_counts_res)
xcell2_counts_res <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/counts/fig1supp_xcell2_counts_res.rds")
benchmark_correlations.counts <- get_xcell2_correlations(benchmark_table = refs.vals.matched.human.xcell2, xCell2results = xcell2_counts_res, round_results = 3, just_xcell2 = TRUE, weight_cors = FALSE, by_val = FALSE)

dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/fpkm"
xcell2_fpkm_res <- paste0(dir2use, "/fig1supp_xcell2_fpkm_res.rds")
get_xcell2_benchmark_results(cyto.vals = human.vals.fpkm, benchmark_table = refs.vals.matched.human.xcell2, save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 6, output_name = xcell2_fpkm_res)
xcell2_fpkm_res <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/fpkm/fig1supp_xcell2_fpkm_res.rds")
benchmark_correlations.fpkm <- get_xcell2_correlations(benchmark_table = refs.vals.matched.human.xcell2, xCell2results = xcell2_fpkm_res, round_results = 3, just_xcell2 = TRUE, weight_cors = FALSE, by_val = FALSE)

benchmark_correlations.tpm$norm <- "TPM"
benchmark_correlations.counts$norm <- "counts"
benchmark_correlations.fpkm$norm <- "FPKM"

benchmark_correlations.allnorm <- rbind(benchmark_correlations.tpm, benchmark_correlations.counts)
benchmark_correlations.allnorm <- rbind(benchmark_correlations.allnorm, benchmark_correlations.fpkm)

saveRDS(benchmark_correlations.allnorm, "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/fpkm/fig1supp_xcell2_all_res.rds")
benchmark_correlations.allnorm <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/fpkm/fig1supp_xcell2_all_res.rds")

fig3a_supp <- benchmark_correlations.allnorm %>%
  ggplot(aes(x = val, y = cor, fill = norm)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Comparison of RNA-Seq Normalization Methods",
    x = "Dataset",
    y = "Spearman rho",
    fill = "Normalization"
  ) +
  # Choose a theme that looks good in print
  theme_minimal(base_size = 14) +
  # Adjust theme elements for clarity
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank()
  )

# fig3a_supp: 1000w x 700h
print(fig3a_supp)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3a_supp.png", plot = fig3a_supp, device = "png", width = 8, height = 6, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3a_supp.pdf", plot = fig3a_supp, device = "pdf", width = 8, height = 6)


# Microarray

human.vals.rma <- human.vals
human.vals.mas5 <- human.vals
human.vals.qn <- human.vals

# , GSE77344 - only oligo
# GSE64385, GSE20300, GSE65135

# GSE64385
GSE64385.rma <- get_expression_matrix(accession = "GSE64385", normalization = "RMA")
GSE64385.mas5 <- get_expression_matrix(accession = "GSE64385", normalization = "MAS5")
GSE64385.qn <- get_expression_matrix(accession = "GSE64385", normalization = "QN")

colnames(GSE64385.rma)
samples <- colnames(human.vals.rma$truth$blood$GSE64385)

colnames(GSE64385.rma) <- gsub("_.*", "", colnames(GSE64385.rma))
colnames(GSE64385.mas5) <- gsub("_.*", "", colnames(GSE64385.mas5))
colnames(GSE64385.qn) <- gsub("_.*", "", colnames(GSE64385.qn))
all(samples %in% colnames(GSE64385.rma))

genes2use <- intersect(intersect(rownames(GSE64385.rma), rownames(GSE64385.mas5)), rownames(GSE64385.qn))
GSE64385.rma <- GSE64385.rma[genes2use,]
GSE64385.mas5 <- GSE64385.mas5[genes2use,]
GSE64385.qn <- GSE64385.qn[genes2use,]
dim(GSE64385.rma)
dim(GSE64385.mas5)
dim(GSE64385.qn)

human.vals.rma$mixtures$blood$GSE64385 <- GSE64385.rma
human.vals.mas5$mixtures$blood$GSE64385 <- GSE64385.mas5
human.vals.qn$mixtures$blood$GSE64385 <- GSE64385.qn

# GSE20300
GSE20300.rma <- get_expression_matrix(accession = "GSE20300", normalization = "RMA")
GSE20300.mas5 <- get_expression_matrix(accession = "GSE20300", normalization = "MAS5")
GSE20300.qn <- get_expression_matrix(accession = "GSE20300", normalization = "QN")

colnames(GSE20300.rma)
samples <- colnames(human.vals.rma$truth$blood$GSE20300)

colnames(GSE20300.rma) <- gsub("\\..*", "", colnames(GSE20300.rma))
colnames(GSE20300.mas5) <- gsub("\\..*", "", colnames(GSE20300.mas5))
colnames(GSE20300.qn) <- gsub("\\..*", "", colnames(GSE20300.qn))
all(samples %in% colnames(GSE20300.rma))

genes2use <- intersect(intersect(rownames(GSE20300.rma), rownames(GSE20300.mas5)), rownames(GSE20300.qn))
GSE20300.rma <- GSE20300.rma[genes2use,]
GSE20300.mas5 <- GSE20300.mas5[genes2use,]
GSE20300.qn <- GSE20300.qn[genes2use,]
dim(GSE20300.rma)
dim(GSE20300.mas5)
dim(GSE20300.qn)

human.vals.rma$mixtures$blood$GSE20300 <- GSE20300.rma
human.vals.mas5$mixtures$blood$GSE20300 <- GSE20300.mas5
human.vals.qn$mixtures$blood$GSE20300 <- GSE20300.qn

# GSE65135
GSE65135.rma <- get_expression_matrix(accession = "GSE65135", normalization = "RMA")
GSE65135.mas5 <- get_expression_matrix(accession = "GSE65135", normalization = "MAS5")
GSE65135.qn <- get_expression_matrix(accession = "GSE65135", normalization = "QN")

colnames(GSE65135.rma)
samples <- colnames(human.vals.rma$truth$blood$GSE65135)

colnames(GSE65135.rma) <- gsub("_.*", "", colnames(GSE65135.rma))
colnames(GSE65135.mas5) <- gsub("_.*", "", colnames(GSE65135.mas5))
colnames(GSE65135.qn) <- gsub("_.*", "", colnames(GSE65135.qn))
all(samples %in% colnames(GSE65135.rma))

genes2use <- intersect(intersect(rownames(GSE65135.rma), rownames(GSE65135.mas5)), rownames(GSE65135.qn))
GSE65135.rma <- GSE65135.rma[genes2use,]
GSE65135.mas5 <- GSE65135.mas5[genes2use,]
GSE65135.qn <- GSE65135.qn[genes2use,]
dim(GSE65135.rma)
dim(GSE65135.mas5)
dim(GSE65135.qn)

human.vals.rma$mixtures$blood$GSE65135 <- GSE65135.rma
human.vals.mas5$mixtures$blood$GSE65135 <- GSE65135.mas5
human.vals.qn$mixtures$blood$GSE65135 <- GSE65135.qn


saveRDS(human.vals.rma, "/bigdata/almogangel/xCell2_data/robust_norm_analysis/human.vals.rma.rds")
saveRDS(human.vals.mas5, "/bigdata/almogangel/xCell2_data/robust_norm_analysis/human.vals.mas5.rds")
saveRDS(human.vals.qn, "/bigdata/almogangel/xCell2_data/robust_norm_analysis/human.vals.qn.rds")


# Run xCell2

ds2check <- c("GSE64385", "GSE20300", "GSE65135")
params2use <- list(
  num_threads = 30,
  min_pb_cells = 20,
  min_pb_samples = 10,
  min_sc_genes = 1e4,
  use_ontology = TRUE,
  return_signatures = FALSE,
  return_analysis = FALSE,
  use_sillover = TRUE,
  spillover_alpha = 0.5
)

refs.vals.matched.human.xcell2 <- refs.vals.matched.human |> 
  mutate(n_val_samples = ncol(human.vals$truth[[val_type]][[val_dataset]])) |> 
  mutate(method = "xCell2", .before = everything()) |> 
  ungroup() |> 
  filter(val_dataset %in% ds2check)


dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/"
xcell2_rma_res <- paste0(dir2use, "/fig1supp_xcell2_rma_res.rds")
get_xcell2_benchmark_results(cyto.vals = human.vals.rma, benchmark_table = refs.vals.matched.human.xcell2, save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 6, output_name = xcell2_rma_res)
xcell2_rma_res <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp//fig1supp_xcell2_rma_res.rds")
benchmark_correlations.rma <- get_xcell2_correlations(benchmark_table = refs.vals.matched.human.xcell2, xCell2results = xcell2_rma_res, round_results = 3, just_xcell2 = TRUE, weight_cors = FALSE, by_val = FALSE)

dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/mas5"
xcell2_mas5_res <- paste0(dir2use, "/fig1supp_xcell2_mas5_res.rds")
get_xcell2_benchmark_results(cyto.vals = human.vals.mas5, benchmark_table = refs.vals.matched.human.xcell2, save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 6, output_name = xcell2_mas5_res)
xcell2_mas5_res <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/mas5/fig1supp_xcell2_mas5_res.rds")
benchmark_correlations.mas5 <- get_xcell2_correlations(benchmark_table = refs.vals.matched.human.xcell2, xCell2results = xcell2_mas5_res, round_results = 3, just_xcell2 = TRUE, weight_cors = FALSE, by_val = FALSE)

dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/qn"
xcell2_qn_res <- paste0(dir2use, "/fig1supp_xcell2_qn_res.rds")
get_xcell2_benchmark_results(cyto.vals = human.vals.qn, benchmark_table = refs.vals.matched.human.xcell2, save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 6, output_name = xcell2_qn_res)
xcell2_qn_res <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/qn/fig1supp_xcell2_qn_res.rds")
benchmark_correlations.qn <- get_xcell2_correlations(benchmark_table = refs.vals.matched.human.xcell2, xCell2results = xcell2_qn_res, round_results = 3, just_xcell2 = TRUE, weight_cors = FALSE, by_val = FALSE)

benchmark_correlations.rma$norm <- "RMA"
benchmark_correlations.mas5$norm <- "MAS5"
benchmark_correlations.qn$norm <- "QN"

benchmark_correlations.allnorm <- rbind(benchmark_correlations.rma, benchmark_correlations.mas5)
benchmark_correlations.allnorm <- rbind(benchmark_correlations.allnorm, benchmark_correlations.qn)

saveRDS(benchmark_correlations.allnorm, "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/fpkm/fig1supp_xcell2_all_res_microarray.rds")
benchmark_correlations.allnorm <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1supp/fpkm/fig1supp_xcell2_all_res_microarray.rds")


fig3b_supp <- benchmark_correlations.allnorm %>%
  ggplot(aes(x = val, y = cor, fill = norm)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Comparison of Microarray Normalization Methods",
    x = "Dataset",
    y = "Spearman rho",
    fill = "Normalization"
  ) +
  # Choose a theme that looks good in print
  theme_minimal(base_size = 14) +
  # Adjust theme elements for clarity
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank()
  )

# fig3a_supp: 1000w x 700h
print(fig3b_supp)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_supp.png", plot = fig3b_supp, device = "png", width = 8, height = 6, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_supp.pdf", plot = fig3b_supp, device = "pdf", width = 8, height = 6)


# Reference affect accuracy -------------------


mouse.refs <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/mouse_references.rds")
mouse.vals <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/mouse_validation.rds")

human.refs <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/human_references.rds")
human.vals <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/human_validation.rds")


params <- list(
  num_threads = 30,
  min_pb_cells = 20,
  min_pb_samples = 10,
  min_sc_genes = 1e4,
  use_ontology = TRUE,
  return_signatures = FALSE,
  return_analysis = FALSE,
  use_sillover = TRUE,
  spillover_alpha = 0.5
)

parallel_param <- MulticoreParam(workers = params$num_threads)


# Human blood vs. Human tumor and Mouse 

ct2use <- "CD8-positive, alpha-beta T cell"
human.vals$truth$blood

refs <- names(human.vals$truth$blood)






