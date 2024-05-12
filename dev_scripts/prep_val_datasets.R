library(tidyverse)
library(limma)
library(GEOquery)
library(affy)
options(timeout=1000)

downloadGEO <- function(geoSeriesId) {

  # get supplementary files
  getGEOSuppFiles(geoSeriesId, baseDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data")

  file_path <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/", geoSeriesId, "/", geoSeriesId, "_RAW.tar")
  untar_path <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/", geoSeriesId, "/data")
  # untar files
  untar(file_path, exdir = untar_path)

  # reading in .cel files
  raw.data <- ReadAffy(celfile.path = untar_path)


  # get expression estimates
  raw.expr <- rma(raw.data, background = FALSE, normalize = FALSE)
  raw.expr <- as.data.frame(exprs(raw.expr))
  raw.expr <- 2^raw.expr

  return(raw.expr)
}

# Prepare validation datasets

celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

load("/bigdata/almogangel/xCell2_data/gene_cov.rda")

# SDY67 (RNA-Seq) ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/SDY67/Bulk_SDY67.RData")
exp.in <- Bulk_SDY67
range(exp.in)

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/SDY67/fSDY67.RData")
frac.in <- fSDY67
frac.in <- frac.in[,-1]

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

# Convert to TPM
exp.in <- GeoTcgaData::countToTpm(exp.in, keyType = "SYMBOL", gene_cov = gene_cov)

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/SDY67_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/SDY67.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# SYD311 ---------------

# For raw data:
# sdy311 <- readRDS("/bigdata/almogangel/microArray_data/sdy311.rds")
# # sdy311.exp <- sdy311$expr
# raw.data <- read.csv("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/SDY311/SDY311_EXP13635_microarray.703343.tsv", sep = "\t")
# metadata <- read.csv("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/SDY311/SDY311-DR50_Subject_2_Illumina_BeadArray.txt", sep = "\t")
# raw.data <- raw.data[!is.na(raw.data$SYMBOL),]
# rownames(raw.data) <- make.unique(raw.data$SYMBOL)
# raw.data <- raw.data[,which(endsWith(x = colnames(raw.data), suffix = "_SIGNAL"))]
# colnames(raw.data) <- gsub("_SIGNAL", "", colnames(raw.data))
# sdy311.fcs <- sdy311$fcs
# metadata <- metadata[metadata$Subject.Accession %in% colnames(sdy311.fcs),c("Subject.Accession", "Expsample.Accession")]
# rownames(metadata) <- metadata$Expsample.Accession
# # Remove duplicates
# dups <- metadata$Subject.Accession[which(duplicated(metadata$Subject.Accession))]
# samples2use <- rownames(metadata[!metadata$Subject.Accession %in% dups,])
# metadata <- metadata[samples2use,]
# raw.data <- raw.data[,samples2use]
# colnames(raw.data) <- metadata[colnames(raw.data),1]
# samples <- intersect(colnames(sdy311.fcs), colnames(raw.data))
# # Remove outlier samples - https://github.com/dviraran/xCell/blob/master/vignettes/xCell-Immport.Rmd
# samples <- samples[!samples %in% c("SUB134240","SUB134283")]
# sdy311.exp <- raw.data
# sdy311.exp <- sdy311.exp[,samples]
# sdy311.fcs <- sdy311.fcs[,samples]
# all(colnames(sdy311.fcs) ==  colnames(sdy311.exp))
# sdy311.fcs <- sdy311.fcs*100
# as.tibble(sdy311.exp[,1:15]) %>% pivot_longer(cols=colnames(sdy311.exp[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()
#
#
# sdy311.exp <- cbind("Gene" = rownames(sdy311.exp), sdy311.exp)
# write.table(sdy311.exp, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/SDY311_expressions.tsv",
#             row.names = FALSE, quote = FALSE, sep = "\t")
#
# rownames(sdy311.fcs) <- plyr::mapvalues(rownames(sdy311.fcs), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
#
# write.table(sdy311.fcs, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/SDY311.tsv",
#             row.names = TRUE, quote = FALSE, sep = "\t")

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/SDY311/Bulk_SDY311.RData")
exp.in <- Bulk_SDY311
# sdy311 <- readRDS("/bigdata/almogangel/microArray_data/sdy311.rds")
# exp.in <- sdy311$expr


range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()
# RMA normalize ?

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/SDY311/fSDY311.RData")
frac.in <- t(fSDY311)
frac.in2 <- sdy311$fcs

samples <- intersect(colnames(frac.in), colnames(exp.in))
colSums(frac.in)

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

frac.in <- frac.in*100
rownames(frac.in) <- rownames(frac.in2)


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/SDY311_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/SDY311.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")


# SYD420  ---------------


sdy420 <- readRDS("/bigdata/almogangel/microArray_data/sdy420.rds")
sdy420.exp <- sdy420$expr
range(sdy420.exp)


sdy420.fcs <- sdy420$fcs
samples <- intersect(colnames(sdy420.exp), colnames(sdy420.fcs))
sdy420.exp <- sdy420.exp[,samples]
sdy420.fcs <- sdy420.fcs[,samples]
all(colnames(sdy420.fcs) ==  colnames(sdy420.exp))

sdy420.fcs <- sdy420.fcs*100


sdy420.exp <- cbind("Gene" = rownames(sdy420.exp), sdy420.exp)
write.table(sdy420.exp, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/SDY420_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(sdy420.fcs) <- plyr::mapvalues(rownames(sdy420.fcs), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(sdy420.fcs, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/SDY420.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE64385  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE64385/Bulk_GSE64385.Rdata")
exp.in <- Bulk_GSE64385
range(exp.in)
colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))
# RMA

# exp.in <- downloadGEO("GSE64385")
# # map probe IDs to gene symbols
# gse <- getGEO("GSE64385", GSEMatrix = TRUE)
# feature.data <- gse[[1]]@featureData@data
# feature.data <- feature.data[,c("ID", "Gene Symbol")]
#
# exp.in <- exp.in %>%
#   rownames_to_column(var = 'ID') %>%
#   inner_join(., feature.data, by = 'ID') %>%
#   separate(`Gene Symbol`, into = c("symbol1", "symbol2"), sep = " /// ") %>%
#   dplyr::select(-ID)
#
# gene.names <- make.unique(exp.in$symbol1)
# exp.in <- as.data.frame(exp.in[,-c(ncol(exp.in), ncol(exp.in)-1)])
# rownames(exp.in) <- gene.names


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE64385/fGSE64385.RData")
frac.in <- fGSE64385
frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names
frac.in <- apply(frac.in, 2, function(x){x/sum(x)})

samples <- intersect(colnames(frac.in), colnames(exp.in))
samples <- samples[!samples %in% c("GSM1570043", "GSM1570044")] # Those sample contain only cancer cells

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

frac.in <- frac.in*100


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE64385_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE64385.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE65133 (*use lumi - fix and run again)  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE65133/Bulk_GSE65133.Rdata")
exp.in <- Bulk_GSE65133
range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()

exp.in <- lumi::lumi.N(exp.in)
# RMA normalization
exp.in <-  limma::backgroundCorrect(as.data.frame(exp.in), method='normexp')
exp.in <- limma::normalizeBetweenArrays(exp.in)
exp.in <- log2(exp.in+1)


# # get supplementary files
# getGEOSuppFiles("GSE65133", baseDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data")
#
# # untar files
# gunzip("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE65133/GSE65133_non-normalized.txt.gz", overwrite = TRUE)
#
# raw.data <- read.csv("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE65133/GSE65133_non-normalized.txt", sep = "\t")
# raw.data <- raw.data[-c(1:2),]
# probes_id <- raw.data[-c(1:2),1]
# raw.data <- raw.data[,which(startsWith(x = as.character(raw.data[1,]), prefix = "GSM"))]
# samples_id <- as.character(raw.data[1,])
# raw.data <- raw.data[-c(1:2),]
# # raw.data <- sapply(raw.data, as.numeric)
# colnames(raw.data) <- samples_id
# rownames(raw.data) <- probes_id
#
# # map probe IDs to gene symbols
# gse <- getGEO("GSE65133", GSEMatrix = TRUE)
# feature.data <- as.data.frame(gse[[1]]@featureData@data)
# feature.data <- feature.data[,c("ID", "Symbol")]
#
# exp.in <- raw.data %>%
#   rownames_to_column(var = 'ID') %>%
#   inner_join(., feature.data, by = 'ID') %>%
#   dplyr::select(-ID)
#
# gene.names <- make.unique(exp.in$Symbol)
# exp.in <- as.data.frame(exp.in[,-ncol(exp.in)])
# rownames(exp.in) <- gene.names



load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE65133/fGSE65133.RData")
frac.in <- fGSE65133

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

frac.in <- frac.in*100


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE65133_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
memcd4 <- colSums(frac.in[5:6,])
frac.in <- frac.in[-c(5:6),]
frac.in <- rbind("CD4-positive, alpha-beta memory T cell" = memcd4, frac.in)
colSums(frac.in)

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE65133.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE106898 (*use lumi - fix and run again)  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE106898/Bulk_GSE106898.RData")
exp.in <- Bulk_GSE106898
range(exp.in)
as.tibble(exp.in[,1:10]) %>% pivot_longer(cols=colnames(exp.in[,1:10])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()

# RMA normalization
exp.in <-  limma::backgroundCorrect(as.data.frame(exp.in), method='normexp')
exp.in <- limma::normalizeBetweenArrays(exp.in)
exp.in <- log2(exp.in+1)
range(exp.in)


colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))
colnames(exp.in) <- gsub("^X", "", colnames(exp.in))

# # get supplementary files
# getGEOSuppFiles("GSE106898", baseDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data")
#
# # untar files
# gunzip("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE106898/GSE106898_non-normalized.txt.gz", overwrite = TRUE)
#
# raw.data <- read.csv("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE106898/GSE106898_non-normalized.txt", sep = "\t")
# probes_id <- raw.data[,1]
# raw.data <- raw.data[,which(endsWith(x = colnames(raw.data), suffix = "_PBMC"))]
# rownames(raw.data) <- probes_id
#
# # map probe IDs to gene symbols
# gse <- getGEO("GSE106898", GSEMatrix = TRUE)
# feature.data <- as.data.frame(gse[[1]]@featureData@data)
# feature.data <- feature.data[,c("ID", "Symbol")]
#
# exp.in <- raw.data %>%
#   rownames_to_column(var = 'ID') %>%
#   inner_join(., feature.data, by = 'ID') %>%
#   dplyr::select(-ID)
#
# gene.names <- make.unique(exp.in$Symbol)
# exp.in <- as.data.frame(exp.in[,-ncol(exp.in)])
# rownames(exp.in) <- gene.names
# colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))
# colnames(exp.in) <- gsub("^X", "", colnames(exp.in))


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE106898/fGSE106898.RData")
frac.in <- fGSE106898

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE106898_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE106898.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE64655 (*RNA-Seq RPKM - fix and run again)  ---------------


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE64655/Bulk_GSE64655.RData")
exp.in <- Bulk_GSE64655
range(exp.in)
exp.in <- 2^exp.in
exp.in <- exp.in-1
colSums(exp.in)

#  Convert to TPM
exp.in <- (exp.in/colSums(exp.in))*10^6
exp.in <- log2(exp.in+1)

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE64655/fGSE64655.RData")
frac.in <- fGSE64655

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE64655_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
frac.in <- frac.in/100

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE64655.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE59654 (*use lumi - fix and run again) ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE59654/Bulk_GSE59654.Rdata")
exp.in <- Bulk_GSE59654
range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()

# RMA normalization
exp.in <-  limma::backgroundCorrect(as.data.frame(exp.in), method='normexp')
exp.in <- limma::normalizeBetweenArrays(exp.in)
exp.in <- log2(exp.in+1)
range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


# colnames(exp.in)

# exp.in <- downloadGEO("GSE59654")
#
# # get supplementary files
# getGEOSuppFiles("GSE59654", baseDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data")
#
# # untar files
# gunzip("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE59654/GSE59654_PBMC.raw.corrected.txt.gz", overwrite = TRUE)
#
# raw.data <- read.csv("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE59654/GSE59654_PBMC.raw.corrected.txt", sep = "\t")
# View(raw.data)
# probes_id <- raw.data[,1]
# raw.data <- raw.data[,grep("PBMC_\\d+$", names(raw.data))]
# rownames(raw.data) <- probes_id
#
# # map probe IDs to gene symbols
# gse <- getGEO("GSE59654", GSEMatrix = TRUE)
# feature.data <- as.data.frame(gse[[1]]@featureData@data)
# feature.data <- feature.data[,c("ID", "Symbol")]
#
# exp.in <- raw.data %>%
#   rownames_to_column(var = 'ID') %>%
#   inner_join(., feature.data, by = 'ID') %>%
#   dplyr::select(-ID)
#
# gene.names <- make.unique(exp.in$Symbol)
# exp.in <- as.data.frame(exp.in[,-ncol(exp.in)])
# rownames(exp.in) <- gene.names


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE59654/fGSE59654.Rdata")
frac.in <- fGSE59654

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names
colnames(exp.in) <- gsub("^X", "", colnames(exp.in))

samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))



exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE59654_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE59654.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE107990  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE107990/Bulk_GSE107990.Rdata")
exp.in <- Bulk_GSE107990
range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


#exp.in <- downloadGEO("GSE107990")

# get supplementary files
# getGEOSuppFiles("GSE107990", baseDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data")
#
# # untar files
# gunzip("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE107990/GSE107990_non-normalized.txt.gz", overwrite = TRUE)
#
# raw.data <- read.csv("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE107990/GSE107990_non-normalized.txt", sep = "\t")
# View(raw.data)
# probes_id <- raw.data[,1]
# raw.data <- raw.data[,grep("AVG_Signal", names(raw.data))]
# rownames(raw.data) <- probes_id
#
#
# # map probe IDs to gene symbols
# gse <- getGEO(filename= "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE107990/GSE107990_series_matrix.txt")
# feature.data <- as.data.frame(gse@featureData@data)
# feature.data <- feature.data[,c("ID", "Symbol")]
#
# exp.in <- raw.data %>%
#   rownames_to_column(var = 'ID') %>%
#   inner_join(., feature.data, by = 'ID') %>%
#   dplyr::select(-ID)
#
# gene.names <- make.unique(exp.in$Symbol)
# exp.in <- as.data.frame(exp.in[,-ncol(exp.in)])
# rownames(exp.in) <- gene.names
#
#
# # /bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE107990/GSE107990.R
# samples2use <- grep("_D0_", colnames(exp.in))
# exp.in <- exp.in[,samples2use]
# sample_names <- sapply(strsplit(colnames(exp.in),"_batch"), "[",1)
# sample_names <- gsub("SNF","X",sample_names)
# sample_names <- gsub("X00","X",sample_names)
# sample_names <- gsub("X0","X",sample_names)
# colnames(exp.in) <- sample_names


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE107990/fGSE107990.Rdata")
frac.in <- fGSE107990

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names

samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

# Remove NAs
# exp.in <- na.omit(exp.in)

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE107990_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE107990.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE20300  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE20300/Bulk_GSE20300.Rdata")
exp.in <- Bulk_GSE20300
range(exp.in)


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE20300/fGSE20300.RData")
frac.in <- fGSE20300

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))




exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE20300_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE20300.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE65135  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE65135/Bulk_GSE65135.Rdata")
exp.in <- Bulk_GSE65135
range(exp.in)

# exp.in <- downloadGEO("GSE65135")
#
# # map probe IDs to gene symbols
# gse <- getGEO("GSE65135", GSEMatrix = TRUE)
# feature.data <- gse[[1]]@featureData@data
# feature.data <- feature.data[,c("ID", "Gene Symbol")]
#
# exp.in <- exp.in %>%
#   rownames_to_column(var = 'ID') %>%
#   inner_join(., feature.data, by = 'ID') %>%
#   separate(`Gene Symbol`, into = c("symbol1", "symbol2"), sep = " /// ") %>%
#   dplyr::select(-ID)
#
# gene.names <- make.unique(exp.in$symbol1)
# exp.in <- as.data.frame(exp.in[,-c(ncol(exp.in), ncol(exp.in)-1)])
# rownames(exp.in) <- gene.names
# colnames(exp.in) <- gsub(".CEL.gz", "", colnames(exp.in))
# colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE65135/fGSE65135.RData")
frac.in <- fGSE65135

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

frac.in <- frac.in*100

as.tibble(exp.in[,1:10]) %>% pivot_longer(cols=colnames(exp.in[,1:10])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE65135_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE65135.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE77343  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE77343/Bulk_GSE77343.RData")
exp.in <- Bulk_GSE77343
range(exp.in)

# geoSeriesId <- "GSE77343"
# file_path <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/", geoSeriesId, "/", geoSeriesId, "_RAW.tar")
# untar_path <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/", geoSeriesId, "/data")
# # untar files
# untar(file_path, exdir = untar_path)
#
# geneCELs <- list.files(untar_path, full = TRUE)
# affyData <- read.celfiles(geneCELs)
# data.eset <- oligo::rma(affyData, background = TRUE, normalize = FALSE)
# exp.in <- exprs(data.eset)
# exp.in <- 2^exp.in
#
#
# # map probe IDs to gene symbols
# GPL <- data.table::fread("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE77344/GPL11532-32230.txt")
# GPL11532 <- GPL[,c(1,10)]
# GPL11532 <- GPL11532[which(GPL11532$gene_assignment!="---"),]
# tmp <- strsplit(GPL11532$gene_assignment,"//")
# gene <- lapply(tmp ,"[",2 )
# gene <- unlist(gene )
# gene <- gsub(" ","",gene)
# GPL11532$gene_assignment <- gene
# symbolID <- sapply(1: length(rownames(exp.in)), function(i){
#   as.character(GPL11532[which(as.character(GPL11532$ID)==rownames(exp.in)[i]),"gene_assignment"])
# })
# symbolID[symbolID == "character(0)"] <- NA
# expr.df <- cbind(symbolID,as.data.frame(exp.in))
# expr.df <- expr.df[!is.na(expr.df$symbolID),]
#
# # Combine duplicated genes
# genes <- expr.df$symbolID
# expr.df <- aggregate(as.matrix(expr.df[,-1])~genes,FUN=mean)
# genes <- expr.df$genes
# expr.df <- expr.df[,-1]
# rownames(expr.df) <- genes
#
# exp.in <- expr.df
# colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE77343/fGSE77343.RData")
frac.in <- fGSE77343

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

frac.in <- frac.in*100


as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE77343_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE77343.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE77344  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE77344/Bulk_GSE77344.RData")
exp.in <- Bulk_GSE77344
range(exp.in)

# geoSeriesId <- "GSE77344"
# file_path <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/", geoSeriesId, "/", geoSeriesId, "_RAW.tar")
# untar_path <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/", geoSeriesId, "/data")
# # untar files
# untar(file_path, exdir = untar_path)
#
# geneCELs <- list.files(untar_path, full = TRUE)
# affyData <- read.celfiles(geneCELs)
# data.eset <- oligo::rma(affyData, background = TRUE, normalize = FALSE)
# exp.in <- exprs(data.eset)
# exp.in <- 2^exp.in
#
#
# # map probe IDs to gene symbols
# GPL <- data.table::fread("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/array_data/GSE77344/GPL11532-32230.txt")
# GPL11532 <- GPL[,c(1,10)]
# GPL11532 <- GPL11532[which(GPL11532$gene_assignment!="---"),]
# tmp <- strsplit(GPL11532$gene_assignment,"//")
# gene <- lapply(tmp ,"[",2 )
# gene <- unlist(gene )
# gene <- gsub(" ","",gene)
# GPL11532$gene_assignment <- gene
# symbolID <- sapply(1: length(rownames(exp.in)), function(i){
#   as.character(GPL11532[which(as.character(GPL11532$ID)==rownames(exp.in)[i]),"gene_assignment"])
# })
# symbolID[symbolID == "character(0)"] <- NA
# expr.df <- cbind(symbolID,as.data.frame(exp.in))
# expr.df <- expr.df[!is.na(expr.df$symbolID),]
#
# # Combine duplicated genes
# genes <- expr.df$symbolID
# expr.df <- aggregate(as.matrix(expr.df[,-1])~genes,FUN=mean)
# genes <- expr.df$genes
# expr.df <- expr.df[,-1]
# rownames(expr.df) <- genes
#
#
# exp.in <- expr.df
# colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))


load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE77344/fGSE77344.RData")
frac.in <- fGSE77344

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

frac.in <- frac.in*100


as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE77344_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE77344.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE93722  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE93722/Bulk_GSE93722.RData")
exp.in <- Bulk_GSE93722
range(exp.in)

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE93722/fGSE93722.RData")
frac.in <- fGSE93722

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE93722_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE93722.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE130824 (RNA-Seq) ---------------

exp.in <- read.table("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/GSE130824/GSE130824_dataNormedFiltered.txt", header = T, row.names = 1, sep = "\t")
range(exp.in)

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/GSE130824/majorCellTypes.rda")
frac.in <- majorCellTypes

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))


exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE130824_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE130824.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

