## ----eval = FALSE-------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
#  }
#  BiocManager::install("xCell2")

## ----eval = FALSE-------------------------------------------------------------
#  if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
#  }
#  devtools::install_github("AlmogAngel/xCell2")

## ----eval = TRUE--------------------------------------------------------------
library(xCell2)

## ----eval = TRUE--------------------------------------------------------------
# Load the demo data
data(dice_demo_ref, package = "xCell2")

# Extract reference matrix
dice_ref <- as.matrix(dice_demo_ref@assays@data$logcounts)
colnames(dice_ref) <- make.unique(colnames(dice_ref))

# Prepare labels data frame
dice_labels <- as.data.frame(dice_demo_ref@colData)
dice_labels$ont <- NA
dice_labels$sample <- colnames(dice_ref)
dice_labels$dataset <- "DICE"

## ----eval = TRUE--------------------------------------------------------------
dice_labels[dice_labels$label == "B cells", ]$ont <- "CL:0000236"
dice_labels[dice_labels$label == "Monocytes", ]$ont <- "CL:0000576"
dice_labels[dice_labels$label == "NK cells", ]$ont <- "CL:0000623"
dice_labels[dice_labels$label == "T cells, CD8+", ]$ont <- "CL:0000625"
dice_labels[dice_labels$label == "T cells, CD4+", ]$ont <- "CL:0000624"
dice_labels[dice_labels$label == "T cells, CD4+, memory", ]$ont <- "CL:0000897"

## ----eval = TRUE--------------------------------------------------------------
xCell2::xCell2GetLineage(labels = dice_labels, outFile = "demo_dice_dep.tsv")

## ----eval = TRUE--------------------------------------------------------------
set.seed(123) # (optional) For reproducibility
DICE_demo.xCell2Ref <- xCell2::xCell2Train(
 ref = dice_ref,
 labels = dice_labels,
 refType = "rnaseq"
)

## ----eval = FALSE-------------------------------------------------------------
#  save(DICE_demo.xCell2Ref, file = "DICE_demo.xCell2Ref.rda")

## ----eval = TRUE--------------------------------------------------------------
data(BlueprintEncode.xCell2Ref)

## ----eval = FALSE-------------------------------------------------------------
#  # Set the URL of the pre-trained reference
#  ref_url <- "https://dviraran.github.io/xCell2refs/references/BlueprintEncode.xCell2Ref.rds"
#  # Set the local filename to save the reference
#  local_filename <- "BlueprintEncode.xCell2Ref.rds"
#  # Download the file
#  download.file(ref_url, local_filename, mode = "wb")
#  # Load the downloaded reference
#  BlueprintEncode.xCell2Ref <- readRDS(local_filename)

## -----------------------------------------------------------------------------
# Load the demo reference object
data(DICE_demo.xCell2Ref, package = "xCell2")
# Load a sample bulk expression dataset
data(mix_demo, package = "xCell2")

## -----------------------------------------------------------------------------
xcell2_results <- xCell2::xCell2Analysis(
 mix = mix_demo,
 xcell2object = DICE_demo.xCell2Ref
)

## -----------------------------------------------------------------------------
sessionInfo()

