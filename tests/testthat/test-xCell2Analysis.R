library(testthat)
library(xCell2)

test_that("xCell2Analysis handles valid input correctly", {
  # Load demo data
  data(DICE_demo.xCell2Ref, package = "xCell2")
  data(mix_demo, package = "xCell2")
  
  # Run analysis
  result <- xCell2Analysis(mix = mix_demo, xcell2object = DICE_demo.xCell2Ref, spillover = FALSE)
  
  # Check output type
  expect_true(inherits(result, "matrix"))
  
  # Check dimensions of the result
  expect_equal(nrow(result), length(unique(gsub("#.*", "", names(getSignatures(DICE_demo.xCell2Ref))))))
  expect_equal(ncol(result), ncol(mix_demo))
  
  # Ensure no negative values
  expect_true(all(result >= 0))
})

test_that("xCell2Analysis handles spillover correction", {
  # Load demo data
  data(DICE_demo.xCell2Ref, package = "xCell2")
  data(mix_demo, package = "xCell2")
  
  # Run analysis with spillover correction
  result <- xCell2Analysis(mix = mix_demo, xcell2object = DICE_demo.xCell2Ref, spillover = TRUE)
  
  # Check output type
  expect_true(inherits(result, "matrix"))
  
  # Ensure no negative values after spillover correction
  expect_true(all(result >= 0))
})

test_that("xCell2Analysis handles insufficient shared genes", {
  # Load demo data
  data(DICE_demo.xCell2Ref, package = "xCell2")
  data(mix_demo, package = "xCell2")
  
  # Remove most of the genes to simulate low overlap
  reduced_mix <- mix_demo[1:5, ]
  
  # Expect error due to insufficient shared genes
  expect_error(
    xCell2Analysis(mix = reduced_mix, xcell2object = DICE_demo.xCell2Ref),
    "This xCell2 reference shares 0 genes with the mixtures and minSharedGenes = 0.9.\nConsider training a new xCell2 reference or adjusting minSharedGenes."
  )
})

test_that("xCell2Analysis handles rawScores parameter", {
  # Load demo data
  data(DICE_demo.xCell2Ref, package = "xCell2")
  data(mix_demo, package = "xCell2")
  
  # Run analysis with rawScores = TRUE
  result <- xCell2Analysis(mix = mix_demo, xcell2object = DICE_demo.xCell2Ref, rawScores = TRUE, spillover = FALSE)
  
  # Check output type
  expect_true(inherits(result, "matrix"))
  
  # Check dimensions of the result
  expect_equal(nrow(result), length(unique(gsub("#.*", "", names(getSignatures(DICE_demo.xCell2Ref))))))
  expect_equal(ncol(result), ncol(mix_demo))
})

test_that("xCell2Analysis respects BPPARAM", {
  # Load demo data
  data(DICE_demo.xCell2Ref, package = "xCell2")
  data(mix_demo, package = "xCell2")
  
  # Run analysis with parallel processing
  library(BiocParallel)
  parallel_param <- MulticoreParam(workers = 2)
  result <- xCell2Analysis(mix = mix_demo, xcell2object = DICE_demo.xCell2Ref, BPPARAM = parallel_param, spillover = FALSE)
  
  # Check output type
  expect_true(inherits(result, "matrix"))
})