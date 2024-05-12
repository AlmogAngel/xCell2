library(tidyverse)

cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

# RNA-Seq --------------------

vals2use <- c("GSE127813", "GSE107011", "ccRCC_cytof_CD45+", "WU_ccRCC_RCCTC", "BG_blood", "GSE120444", "GSE115823",
              "GSE121127", "NSCLC_cytof", "SDY67", "GSE107572", "GSE53655", "GSE64655", "GSE93722", "GSE130824")
vals.blood <- cyto.vals$mixtures$blood[names(cyto.vals$mixtures$blood) %in% vals2use]
truths.blood <- cyto.vals$truth$blood[names(cyto.vals$truth$blood) %in% vals2use]
vals.tumor <- cyto.vals$mixtures$tumor[names(cyto.vals$mixtures$tumor) %in% vals2use]
truths.tumor <- cyto.vals$truth$tumor[names(cyto.vals$truth$tumor) %in% vals2use]
vals.others <- cyto.vals$mixtures$other[names(cyto.vals$mixtures$other) %in% vals2use]
truths.others <- cyto.vals$truth$other[names(cyto.vals$truth$other) %in% vals2use]

vals <- c(vals.blood, vals.tumor, vals.others)
truths <- c(truths.blood, truths.tumor, truths.others)

# Make sure truth mats are in 0-1 range
for (i in 1:length(truths)) {
  if (max(truths[[i]], na.rm = TRUE) > 1) {
    print(names(truths)[i])
    truths[[i]] <- truths[[i]]/100
  }
}

# # Get shared genes
# genes2use <- Reduce(intersect, lapply(vals, rownames))
# vals <- lapply(vals, function(x){x[genes2use,]})
# for (i in 1:length(vals)) {
#   colnames(vals[[i]]) <- paste0(vals2use[i], "#", colnames(vals[[i]]))
# }

# Combine mixtures
# vals.combined <- do.call(cbind, vals)

# # Get all cell types
# allcts <- unique(unlist(lapply(truths, rownames)))
# for (i in 1:length(truths)) {
#   colnames(truths[[i]]) <- paste0(vals2use[i], "#", colnames(truths[[i]]))
#   cts.missing <- allcts[!allcts %in% rownames(truths[[i]])]
#   tmp_mat <- matrix(data = NA, nrow = length(cts.missing), ncol = ncol(truths[[i]]), dimnames = list(cts.missing, colnames(truths[[i]])))
#   truths[[i]] <- rbind(truths[[i]], tmp_mat)[allcts,]
# }

# # Combine truths
# truths.combined <- do.call(cbind, truths)

# # Check data overlap
# samples <- intersect(colnames(truths.combined), colnames(vals.combined))
# truths.combined <- truths.combined[,samples]
# vals.combined <- vals.combined[,samples]
# all(colnames(truths.combined) == colnames(vals.combined))



val.list <- list(mixture = vals, truth = truths)

saveRDS(val.list, "/bigdata/almogangel/xCell2/data/rnaseq_filtering_data.rds")



# Microarray --------------------

vals2use <- c("SDY311", "SDY420", "GSE64385", "GSE65133", "GSE106898", "GSE59654", "GSE107990", "GSE20300", "GSE65135", "GSE77343", "GSE77344")
vals <- cyto.vals$mixtures$blood[vals2use]
truths <- cyto.vals$truth$blood[vals2use]

# Make sure truth mats are in 0-100% range
for (i in 1:length(truths)) {
  if (max(truths[[i]], na.rm = TRUE) > 1) {
    truths[[i]] <- truths[[i]]/100
  }
}

val.list <- list(mixture = vals, truth = truths)

saveRDS(val.list, "/bigdata/almogangel/xCell2/data/array_filtering_data.rds")





