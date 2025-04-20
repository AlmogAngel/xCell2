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

# For common cell type annotation 
celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_dev/paper_figures/data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

load("/bigdata/almogangel/xCell2_data/gene_cov.rda")

# GSE64655 - can't find (don't use Kassandra's)


# Kassandra RNA-Seq collection (x11 datasets) --------------------

# Those are bulk RNA-Seq TPM normalized validation curated by Zaitsev A (2022):
# BG_blood, GSE107011, GSE107572, GSE115823, GSE120444, GSE127813, GSE64655, NSCLC_cytof, SDY67, WU_ccRCC_RCCTC, ccRCC_cytof_CD45+
# /bigdata/almogangel/kassandra_data/24_validation_datasets/


# SDY67 (RNA-Seq) - Use the one from Kassandra instead ---------------

truth <- as.matrix(read.table("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/SDY67.tsv", header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
rows <- rownames(truth)
truth <- suppressWarnings(apply(truth, 2, as.numeric))
rownames(truth) <- rows
truth <- truth[!apply(truth, 1, function(x){all(is.na(x))}),]
truth <- truth[,!is.na(colSums(truth))]

exp <- as.matrix(read.table("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/SDY67_expressions.tsv", header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
samples2use <- intersect(colnames(truth), colnames(exp))
exp <- exp[,samples2use]
truth <- truth[,samples2use]

exp.in <- cbind("Gene" = rownames(exp), exp)
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/SDY67_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(truth, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/SDY67.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# SYD311 ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/SDY311/Bulk_SDY311.RData")
exp.in <- Bulk_SDY311

range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()
# RMA normalize

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/SDY311/fSDY311.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/SDY311_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/SDY311.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")


# SYD420  ---------------


sdy420 <- readRDS("/bigdata/almogangel/microArray_data/sdy420.rds")
sdy420.exp <- sdy420$expr
range(sdy420.exp)

range(sdy420.exp)
as.tibble(sdy420.exp[,1:15]) %>% pivot_longer(cols=colnames(sdy420.exp[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()
# RMA normalize

sdy420.fcs <- sdy420$fcs
samples <- intersect(colnames(sdy420.exp), colnames(sdy420.fcs))
sdy420.exp <- sdy420.exp[,samples]
sdy420.fcs <- sdy420.fcs[,samples]
all(colnames(sdy420.fcs) ==  colnames(sdy420.exp))

sdy420.fcs <- sdy420.fcs*100


sdy420.exp <- cbind("Gene" = rownames(sdy420.exp), sdy420.exp)
write.table(sdy420.exp, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/SDY420_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(sdy420.fcs) <- plyr::mapvalues(rownames(sdy420.fcs), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(sdy420.fcs, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/SDY420.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE64385  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE64385/Bulk_GSE64385.Rdata")
exp.in <- Bulk_GSE64385
range(exp.in)
colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))

range(exp.in)
as.tibble(exp.in[,1:5]) %>% pivot_longer(cols=colnames(exp.in[,1:5])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()
# RMA normalize


load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE64385/fGSE64385.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE64385_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE64385.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE65133  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE65133/Bulk_GSE65133.Rdata")
exp.in <- Bulk_GSE65133


range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()
# Not RMA


# RMA normalization
exp.in <-  limma::backgroundCorrect(as.data.frame(exp.in), method='normexp')
exp.in <- limma::normalizeBetweenArrays(exp.in)
exp.in <- log2(exp.in+1)

range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()
# Not RMA



load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE65133/fGSE65133.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE65133_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
memcd4 <- colSums(frac.in[5:6,])
frac.in <- frac.in[-c(5:6),]
frac.in <- rbind("CD4-positive, alpha-beta memory T cell" = memcd4, frac.in)
colSums(frac.in)

write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE65133.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE106898 ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE106898/Bulk_GSE106898.RData")
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


load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE106898/fGSE106898.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE106898_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE106898.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE59654  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE59654/Bulk_GSE59654.Rdata")
exp.in <- Bulk_GSE59654
range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()

# RMA normalization
exp.in <-  limma::backgroundCorrect(as.data.frame(exp.in), method='normexp')
exp.in <- limma::normalizeBetweenArrays(exp.in)
exp.in <- log2(exp.in+1)
range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()



load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE59654/fGSE59654.Rdata")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE59654_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE59654.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE107990  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE107990/Bulk_GSE107990.Rdata")
exp.in <- Bulk_GSE107990
range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE107990/fGSE107990.Rdata")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE107990_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE107990.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE20300  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE20300/Bulk_GSE20300.Rdata")
exp.in <- Bulk_GSE20300

range(exp.in)
as.tibble(exp.in[,1:5]) %>% pivot_longer(cols=colnames(exp.in[,1:5])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE20300/fGSE20300.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE20300_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE20300.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE65135  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE65135/Bulk_GSE65135.Rdata")
exp.in <- Bulk_GSE65135

range(exp.in)
as.tibble(exp.in[,1:5]) %>% pivot_longer(cols=colnames(exp.in[,1:5])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE65135/fGSE65135.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE65135_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE65135.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE77343  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE77343/Bulk_GSE77343.RData")
exp.in <- Bulk_GSE77343

range(exp.in)
as.tibble(exp.in[,1:5]) %>% pivot_longer(cols=colnames(exp.in[,1:5])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE77343/fGSE77343.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE77343_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE77343.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE77344  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE77344/Bulk_GSE77344.RData")
exp.in <- Bulk_GSE77344

range(exp.in)
as.tibble(exp.in[,1:15]) %>% pivot_longer(cols=colnames(exp.in[,1:15])) %>% ggplot(.,aes(x=name,y=value)) + geom_boxplot()


load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE77344/fGSE77344.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE77344_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE77344.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE93722 (RNA-Seq) ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE93722/Bulk_GSE93722.RData")
exp.in <- Bulk_GSE93722
range(exp.in)

#  Convert to CPM
exp.in <- apply(exp.in, 2, function(x){x*10^6/sum(x)})
exp.in <- log2(exp.in+1)
colSums(exp.in)

load("/bigdata/almogangel/xCell2_data/benchmarking_data/archive/validations_raw/GSE93722/fGSE93722.RData")
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE93722_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE93722.tsv",
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
write.table(exp.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/GSE130824_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/GSE130824.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")


# petitprez (mouse) -----------
petitprez.exp <- readRDS("/bigdata/almogangel/xCell2_data/mouse/petitprez/petitprez_tpm.rds")
petitprez.exp <- as.matrix(petitprez.exp)

petitprez.exp <- apply(petitprez.exp, 2, function(x){x*10^6/sum(x)})
colSums(petitprez.exp)
petitprez.exp <- log2(petitprez.exp+1)
colSums(petitprez.exp)

petitprez.exp <- cbind("Gene" = rownames(petitprez.exp), petitprez.exp)
write.table(petitprez.exp, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/mouse/petitprez_expressions.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

petitprez.truth <- readRDS("/bigdata/almogangel/xCell2_data/mouse/petitprez/petitprez_facs.rds")
celltype <- rownames(petitprez.truth)
celltype <- plyr::mapvalues(celltype, celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
celltype[celltype  == "MAST cells"] <- "mast cell"
petitprez.truth <- as.matrix(apply(petitprez.truth, 2, as.numeric))
rownames(petitprez.truth) <- celltype
petitprez.truth <- petitprez.truth[!rownames(petitprez.truth) %in% c("Monocytes/Macrophages", "NK/T cells", "B derived cells", "B cells germinal"),]
write.table(petitprez.truth, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/mouse/petitprez.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


# chen (mouse) ------------
chen.exp <- readRDS("/bigdata/almogangel/xCell2_data/mouse/chen/chen_tpm.rds")

range(chen.exp)
chen.exp <- apply(chen.exp, 2, function(x){x*10^6/sum(x)})
colSums(chen.exp)
chen.exp <- log2(chen.exp+1)
colSums(chen.exp)


chen.exp <- cbind("Gene" = rownames(chen.exp), chen.exp)
write.table(chen.exp, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/mouse/chen_expressions.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


chen.truth <- readRDS("/bigdata/almogangel/xCell2_data/mouse/chen/chen_facs.rds")
rownames(chen.truth) <- plyr::mapvalues(rownames(chen.truth), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
rownames(chen.truth)[1] <- "monocyte"
rownames(chen.truth)[2] <- "B cell"
rownames(chen.truth)[4] <- "CD4-positive, alpha-beta T cell"

t_truth <- colSums(chen.truth[c("CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"),])
chen.truth <- rbind(chen.truth, "T cell" = t_truth)

write.table(chen.truth, "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/mouse/chen.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
