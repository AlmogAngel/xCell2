library(testthat)
library(xCell2)

test_that("xCell2Train works with valid RNA-Seq data", {
  # Load demo reference and metadata
  data(dice_demo_ref, package = "xCell2")
  dice_ref <- SummarizedExperiment::assay(dice_demo_ref, "logcounts")
  dice_labels <- as.data.frame(SummarizedExperiment::colData(dice_demo_ref))
  
  # Add ontology and sample columns
  dice_labels$ont <- NA
  dice_labels$sample <- colnames(dice_ref)
  dice_labels$dataset <- "DICE"
  
  # Test RNA-Seq training
  trained_ref <- xCell2Train(ref = dice_ref, labels = dice_labels, refType = "rnaseq")
  
  # Check output type
  expect_s4_class(trained_ref, "xCell2Object")
  
  # Check that all necessary components are present in the S4 object
  expect_true(!is.null(trained_ref@signatures))
  expect_true(!is.null(trained_ref@params))
  expect_true(!is.null(trained_ref@spill_mat))
  expect_true(!is.null(trained_ref@genes_used))
  
  # Ensure there are at least 3 signatures per cell type
  cell_types <- unique(dice_labels$label)
  for (ct in cell_types) {
    sigs <- trained_ref@signatures[startsWith(names(trained_ref@signatures), paste0(ct, "#"))]
    expect_true(length(sigs) >= 3)
  }
})

test_that("xCell2Train works with scRNA-Seq data", {
  # Load demo reference and metadata
  data(dice_demo_ref, package = "xCell2")
  dice_ref <- SummarizedExperiment::assay(dice_demo_ref, "logcounts")
  dice_labels <- as.data.frame(SummarizedExperiment::colData(dice_demo_ref))
  
  # Add ontology and sample columns
  dice_labels$ont <- NA
  dice_labels$sample <- colnames(dice_ref)
  dice_labels$dataset <- "DICE"
  
  # Test scRNA-Seq training
  trained_ref <- xCell2Train(ref = dice_ref, labels = dice_labels, refType = "sc", minScGenes = 800)
  
  # Check output type
  expect_s4_class(trained_ref, "xCell2Object")
  
})

test_that("xCell2Train handles missing ontology gracefully", {
  # Load demo reference and metadata
  data(dice_demo_ref, package = "xCell2")
  dice_ref <- SummarizedExperiment::assay(dice_demo_ref, "logcounts")
  dice_labels <- as.data.frame(SummarizedExperiment::colData(dice_demo_ref))
  
  # Add sample and dataset columns without ontology
  dice_labels$ont <- NA
  dice_labels$sample <- colnames(dice_ref)
  dice_labels$dataset <- "DICE"
  
  # Train reference without ontology
  trained_ref <- xCell2Train(ref = dice_ref, labels = dice_labels, refType = "rnaseq", useOntology = FALSE)
  
  # Check output type
  expect_s4_class(trained_ref, "xCell2Object")
  
  # Ensure dependencies are empty
  expect_true(length(trained_ref@dependencies) == 0)
})

test_that("xCell2Train validates inputs correctly", {
  # Load demo reference and metadata
  data(dice_demo_ref, package = "xCell2")
  dice_ref <- SummarizedExperiment::assay(dice_demo_ref, "logcounts")
  dice_labels <- as.data.frame(SummarizedExperiment::colData(dice_demo_ref))
  
  # Add sample and dataset columns without ontology
  dice_labels$ont <- NA
  dice_labels$sample <- colnames(dice_ref)
  dice_labels$dataset <- "DICE"
  
  # Test with invalid refType
  expect_error(xCell2Train(ref = dice_ref, labels = dice_labels, refType = "invalid"),
               "refType should be 'rnaseq', 'array' or 'sc'")
  
  # Test with insufficient cell types
  insufficient_labels <- dice_labels[1:2, ]
  expect_error(xCell2Train(ref = dice_ref, labels = insufficient_labels, refType = "rnaseq"),
               "Reference must have at least 3 cell types!")
})

test_that("xCell2Train can return analysis results", {
  # Load demo reference and metadata
  data(dice_demo_ref, package = "xCell2")
  dice_ref <- SummarizedExperiment::assay(dice_demo_ref, "logcounts")
  dice_labels <- as.data.frame(SummarizedExperiment::colData(dice_demo_ref))
  
  # Add sample and dataset columns
  dice_labels$ont <- NA
  dice_labels$sample <- colnames(dice_ref)
  dice_labels$dataset <- "DICE"
  
  # Generate a mock bulk mixture
  mock_mix <- dice_ref[, 1:5] + matrix(rnorm(nrow(dice_ref) * 5, sd = 0.1), nrow = nrow(dice_ref))
  
  # Train and return analysis
  analysis_result <- xCell2Train(ref = dice_ref, mix = mock_mix, labels = dice_labels, refType = "rnaseq", returnAnalysis = TRUE)
  
  # Check output type
  expect_true(inherits(analysis_result, "matrix"))
  
  # Ensure dimensions match the mock mixture
  expect_equal(ncol(analysis_result), ncol(mock_mix))
})