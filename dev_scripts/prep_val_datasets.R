library(tidyverse)

# Prepare validation datasets

celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# SDY67  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/SDY67/Bulk_SDY67.RData")
exp.in <- Bulk_SDY67
colnames(exp.in)

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

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/SDY67_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/SDY67.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# SYD311 ---------------
sdy311 <- readRDS("/bigdata/almogangel/microArray_data/sdy311.rds")
sdy311.exp <- sdy311$expr
sdy311.fcs <- sdy311$fcs
samples <- intersect(colnames(sdy311.exp), colnames(sdy311.fcs))
# Remove outlier samples - https://github.com/dviraran/xCell/blob/master/vignettes/xCell-Immport.Rmd
samples <- samples[!samples %in% c("SUB134240","SUB134283")]
sdy311.exp <- sdy311.exp[,samples]
sdy311.fcs <- sdy311.fcs[,samples]
all(colnames(sdy311.fcs) ==  colnames(sdy311.exp))

sdy311.exp <- cbind("Gene" = rownames(sdy311$expr), sdy311.exp)
write.table(sdy311.exp, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/SDY311_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(sdy311.fcs) <- plyr::mapvalues(rownames(sdy311.fcs), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(sdy311.fcs, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/SDY311.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")



# SYD420  ---------------
sdy420 <- readRDS("/bigdata/almogangel/microArray_data/sdy420.rds")
sdy420.exp <- sdy420$expr
sdy420.fcs <- sdy420$fcs
samples <- intersect(colnames(sdy420.exp), colnames(sdy420.fcs))
sdy420.exp <- sdy420.exp[,samples]
sdy420.fcs <- sdy420.fcs[,samples]
all(colnames(sdy420.fcs) ==  colnames(sdy420.exp))

sdy420.exp <- cbind("Gene" = rownames(sdy420.exp), sdy420.exp)
write.table(sdy420.exp, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/SDY420_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(sdy420.fcs) <- plyr::mapvalues(rownames(sdy420.fcs), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)

write.table(sdy420.fcs, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/SDY420.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE64385  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE64385/Bulk_GSE64385.Rdata")
exp.in <- Bulk_GSE64385
colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE64385/fGSE64385.RData")
frac.in <- fGSE64385
frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names
frac.in <- apply(frac.in, 2, function(x){x/sum(x)})


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE64385_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE64385.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE65133  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE65133/Bulk_GSE65133.Rdata")
exp.in <- Bulk_GSE65133
colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))

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

# GSE106898  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE106898/Bulk_GSE106898.RData")
exp.in <- Bulk_GSE106898
colnames(exp.in) <- gsub("_.*", "", colnames(exp.in))
colnames(exp.in) <- gsub("^X", "", colnames(exp.in))

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
frac.in <- frac.in/100

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE106898.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE64655  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE64655/Bulk_GSE64655.RData")
exp.in <- Bulk_GSE64655


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

# GSE59654  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE59654/Bulk_GSE59654.Rdata")
exp.in <- Bulk_GSE59654
colnames(exp.in)

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE59654/fGSE59654.Rdata")
frac.in <- fGSE59654

frac.in <- t(frac.in)
names <- colnames(frac.in)
frac.in <- t(apply(frac.in, 1, as.numeric))
colnames(frac.in) <- names


samples <- intersect(colnames(frac.in), colnames(exp.in))

exp.in <- exp.in[,samples]
frac.in <- frac.in[,samples]
all(colnames(frac.in) ==  colnames(exp.in))

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE59654_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
frac.in <- frac.in/100

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE59654.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE107990  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE107990/Bulk_GSE107990.Rdata")
exp.in <- Bulk_GSE107990
colnames(exp.in)

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

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE107990_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
frac.in <- frac.in/100

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE107990.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE20300  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE20300/Bulk_GSE20300.Rdata")
exp.in <- Bulk_GSE20300
colnames(exp.in)

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
frac.in <- frac.in/100

write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE20300.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE65135  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE65135/Bulk_GSE65135.Rdata")
exp.in <- Bulk_GSE65135
colnames(exp.in)

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

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE65135_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE65135.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE77343  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE77343/Bulk_GSE77343.RData")
exp.in <- Bulk_GSE77343
colnames(exp.in)

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

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE77343_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE77343.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE77344  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE77344/Bulk_GSE77344.RData")
exp.in <- Bulk_GSE77344
colnames(exp.in)

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

exp.in <- cbind("Gene" = rownames(exp.in), exp.in)
write.table(exp.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE77344_expressions.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

rownames(frac.in) <- plyr::mapvalues(rownames(frac.in), celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)


write.table(frac.in, "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE77344.tsv",
            row.names = TRUE, quote = FALSE, sep = "\t")

# GSE93722  ---------------

load("/bigdata/almogangel/xCell2_data/benchmarking_data/validations/new_benchmark_data/data/GSE93722/Bulk_GSE93722.RData")
exp.in <- Bulk_GSE93722
colnames(exp.in)

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
