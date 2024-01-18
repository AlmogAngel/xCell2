library(tidyverse)

cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

# Blood - RNA-Seq - cytometry --------------------

vals2use <- c("GSE127813", "GSE107011", "BG_blood", "SDY67")
vals <- cyto.vals$mixtures$blood[vals2use]
truths <- cyto.vals$truth$blood[vals2use]

# Make sure truth mats are in 0-100% range
for (i in 1:length(truths)) {
  if (max(truths[[i]], na.rm = TRUE) <= 1) {
    truths[[i]] <- truths[[i]]*100
  }
}

# Get shared genes
genes2use <- Reduce(intersect, lapply(vals, rownames))
vals <- lapply(vals, function(x){x[genes2use,]})
for (i in 1:length(vals)) {
  colnames(vals[[i]]) <- paste0(vals2use[i], "#", colnames(vals[[i]]))
}

# Combine mixtures
vals.combined <- do.call(cbind, vals)

# Get all cell types
allcts <- unique(unlist(lapply(truths, rownames)))
for (i in 1:length(truths)) {
  colnames(truths[[i]]) <- paste0(vals2use[i], "#", colnames(truths[[i]]))
  cts.missing <- allcts[!allcts %in% rownames(truths[[i]])]
  tmp_mat <- matrix(data = NA, nrow = length(cts.missing), ncol = ncol(truths[[i]]), dimnames = list(cts.missing, colnames(truths[[i]])))
  truths[[i]] <- rbind(truths[[i]], tmp_mat)[allcts,]
}

# Combine truths
truths.combined <- do.call(cbind, truths)

# Check data overlap
samples <- intersect(colnames(truths.combined), colnames(vals.combined))
truths.combined <- truths.combined[,samples]
vals.combined <- vals.combined[,samples]
all(colnames(truths.combined) == colnames(vals.combined))

# Rank mixture using SingScore to save space
# vals.combined <- singscore::rankGenes(vals.combined)

val.list <- list(mixture = vals.combined, truth = truths.combined)

saveRDS(val.list, "/bigdata/almogangel/xCell2/data/blood_rnaseq_filtering_data.rds")


# Tumor - RNA-Seq - cytometry + others  --------------------

vals2use <- c("ccRCC_cytof_CD45+", "WU_ccRCC_RCCTC", "GSE121127", "NSCLC_cytof")
vals <- cyto.vals$mixtures$tumor[vals2use]
truths <- cyto.vals$truth$tumor[vals2use]

# Make sure truth mats are in 0-100% range
for (i in 1:length(truths)) {
  if (max(truths[[i]], na.rm = TRUE) <= 1) {
    truths[[i]] <- truths[[i]]*100
  }
}

# Get shared genes
genes2use <- Reduce(intersect, lapply(vals, rownames))
vals <- lapply(vals, function(x){x[genes2use,]})
for (i in 1:length(vals)) {
  colnames(vals[[i]]) <- paste0(vals2use[i], "#", colnames(vals[[i]]))
}

# Combine mixtures
vals.combined <- do.call(cbind, vals)

# Get all cell types
allcts <- unique(unlist(lapply(truths, rownames)))
for (i in 1:length(truths)) {
  colnames(truths[[i]]) <- paste0(vals2use[i], "#", colnames(truths[[i]]))
  cts.missing <- allcts[!allcts %in% rownames(truths[[i]])]
  tmp_mat <- matrix(data = NA, nrow = length(cts.missing), ncol = ncol(truths[[i]]), dimnames = list(cts.missing, colnames(truths[[i]])))
  truths[[i]] <- rbind(truths[[i]], tmp_mat)[allcts,]
}

# Combine truths
truths.combined <- do.call(cbind, truths)

# Check data overlap
samples <- intersect(colnames(truths.combined), colnames(vals.combined))
truths.combined <- truths.combined[,samples]
vals.combined <- vals.combined[,samples]
all(colnames(truths.combined) == colnames(vals.combined))

# Rank mixture using SingScore to save space
# vals.combined <- singscore::rankGenes(vals.combined)

val.list <- list(mixture = vals.combined, truth = truths.combined)

saveRDS(val.list, "/bigdata/almogangel/xCell2/data/tumor_rnaseq_filtering_data.rds")

# Blood - Microarray - cytometry --------------------

vals2use <- c("SDY311", "SDY420", "GSE64385", "GSE65133", "GSE106898", "GSE59654", "GSE107990", "GSE65135", "GSE77343", "GSE77344")
vals <- cyto.vals$mixtures$blood[vals2use]
truths <- cyto.vals$truth$blood[vals2use]

# Make sure truth mats are in 0-100% range
for (i in 1:length(truths)) {
  if (max(truths[[i]], na.rm = TRUE) <= 1) {
    truths[[i]] <- truths[[i]]*100
  }
}

# Get shared genes
genes2use <- Reduce(intersect, lapply(vals, rownames))
vals <- lapply(vals, function(x){x[genes2use,]})
for (i in 1:length(vals)) {
  colnames(vals[[i]]) <- paste0(vals2use[i], "#", colnames(vals[[i]]))
}

# Combine mixtures
vals.combined <- do.call(cbind, vals)

# Get all cell types
allcts <- unique(unlist(lapply(truths, rownames)))
for (i in 1:length(truths)) {
  colnames(truths[[i]]) <- paste0(vals2use[i], "#", colnames(truths[[i]]))
  cts.missing <- allcts[!allcts %in% rownames(truths[[i]])]
  tmp_mat <- matrix(data = NA, nrow = length(cts.missing), ncol = ncol(truths[[i]]), dimnames = list(cts.missing, colnames(truths[[i]])))
  truths[[i]] <- rbind(truths[[i]], tmp_mat)[allcts,]
}

# Combine truths
truths.combined <- do.call(cbind, truths)

# Check data overlap
samples <- intersect(colnames(truths.combined), colnames(vals.combined))
truths.combined <- truths.combined[,samples]
vals.combined <- vals.combined[,samples]
all(colnames(truths.combined) == colnames(vals.combined))

# Rank mixture using SingScore to save space
# vals.combined <- singscore::rankGenes(vals.combined)

val.list <- list(mixture = vals.combined, truth = truths.combined)

saveRDS(val.list, "/bigdata/almogangel/xCell2/data/blood_array_filtering_data.rds")





# Blood - combined - cytometry --------------------

vals2use <- c("GSE127813", "GSE107011", "BG_blood", "SDY67",
              "SDY311", "SDY420", "GSE64385", "GSE65133", "GSE106898", "GSE59654", "GSE107990", "GSE65135", "GSE77343", "GSE77344")
vals <- cyto.vals$mixtures$blood[vals2use]
truths <- cyto.vals$truth$blood[vals2use]

# Make sure truth mats are in 0-100% range
for (i in 1:length(truths)) {
  if (max(truths[[i]], na.rm = TRUE) <= 1) {
    truths[[i]] <- truths[[i]]*100
  }
}

# Get shared genes
genes2use <- Reduce(intersect, lapply(vals, rownames))
vals <- lapply(vals, function(x){x[genes2use,]})
for (i in 1:length(vals)) {
  colnames(vals[[i]]) <- paste0(vals2use[i], "#", colnames(vals[[i]]))
}

# Combine mixtures
vals.combined <- do.call(cbind, vals)

# Get all cell types
allcts <- unique(unlist(lapply(truths, rownames)))
for (i in 1:length(truths)) {
  colnames(truths[[i]]) <- paste0(vals2use[i], "#", colnames(truths[[i]]))
  cts.missing <- allcts[!allcts %in% rownames(truths[[i]])]
  tmp_mat <- matrix(data = NA, nrow = length(cts.missing), ncol = ncol(truths[[i]]), dimnames = list(cts.missing, colnames(truths[[i]])))
  truths[[i]] <- rbind(truths[[i]], tmp_mat)[allcts,]
}

# Combine truths
truths.combined <- do.call(cbind, truths)

# Check data overlap
samples <- intersect(colnames(truths.combined), colnames(vals.combined))
truths.combined <- truths.combined[,samples]
vals.combined <- vals.combined[,samples]
all(colnames(truths.combined) == colnames(vals.combined))

# Rank mixture using SingScore to save space
# vals.combined <- singscore::rankGenes(vals.combined)

val.list <- list(mixture = vals.combined, truth = truths.combined)

saveRDS(val.list, "/bigdata/almogangel/xCell2/data/blood_filtering_data.rds")


