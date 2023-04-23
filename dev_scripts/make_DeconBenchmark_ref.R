library(DeconBenchmark)

# Load blood super ref labels
labels <- readRDS("/bigdata/almogangel/xCell2/dev_data/sref_blood_labels.rds")
ref <- readRDS("/bigdata/almogangel/xCell2/dev_data/sref_blood_data.rds")

singleCellExpr <- ref
singleCellLabels <- labels$label
singleCellSubjects <- labels$dataset

# Generate reference
reference <- generateReference(singleCellExpr, singleCellLabels, types = c("markers", "sigGenes", "signature", "cellTypeExpr"))

# Save reference
saveRDS(reference, "/bigdata/almogangel/xCell2/dev_data/DeconBenchmark_sref_blood_reference.rds")
