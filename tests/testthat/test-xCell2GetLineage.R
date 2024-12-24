library(testthat)
library(xCell2)

test_that("xCell2GetLineage handles valid input correctly", {
  # Mock labels input
  labels <- data.frame(
    ont = c("CL:0000236", "CL:0000576", "CL:0000623", NA, NA),
    label = c("B cells", "Monocytes", "NK cells", "Unknown1", "Unknown2"),
    sample = c("S1", "S2", "S3", "S4", "S5"),
    dataset = "Demo"
  )
  
  # Call the function without output file
  lineage <- xCell2GetLineage(labels)
  
  # Check that the output is a list with the correct names
  expect_type(lineage, "list")
  expect_named(lineage, c("B cells", "Monocytes", "NK cells", "Unknown1", "Unknown2"))
  
  # Check that descendants and ancestors are character vectors
  expect_type(lineage[["B cells"]]$descendants, "character")
  expect_type(lineage[["B cells"]]$ancestors, "character")
})

test_that("xCell2GetLineage handles missing ontologies gracefully", {
  # Mock labels input with no ontology information
  labels <- data.frame(
    ont = c(NA, NA, NA),
    label = c("Cell1", "Cell2", "Cell3"),
    sample = c("S1", "S2", "S3"),
    dataset = "Demo"
  )
  
  # Call the function
  lineage <- xCell2GetLineage(labels)
  
  # Check that the output is a list with no dependencies
  expect_type(lineage, "list")
  expect_named(lineage, c("Cell1", "Cell2", "Cell3"))
  expect_equal(lineage[["Cell1"]]$descendants, character(0))
  expect_equal(lineage[["Cell1"]]$ancestors, character(0))
})

test_that("xCell2GetLineage handles file output correctly", {
  # Mock labels input
  labels <- data.frame(
    ont = c("CL:0000236", "CL:0000576", NA),
    label = c("B cells", "Monocytes", "Unknown"),
    sample = c("S1", "S2", "S3"),
    dataset = "Demo"
  )
  
  # Specify an output file
  temp_file <- tempfile(fileext = ".tsv")
  
  # Call the function with file output
  xCell2GetLineage(labels, outFile = temp_file)
  
  # Check that the file is created
  expect_true(file.exists(temp_file))
  
  # Check the contents of the file
  content <- readr::read_tsv(temp_file)
  expect_equal(ncol(content), 4) # ont, label, descendants, ancestors
  expect_equal(nrow(content), 3) # Same as number of input labels
  
  # Clean up
  unlink(temp_file)
})

test_that("xCell2GetLineage handles empty input gracefully", {
  # Empty labels input
  labels <- data.frame(
    ont = character(0),
    label = character(0),
    sample = character(0),
    dataset = character(0)
  )
  
  # Call the function
  lineage <- xCell2GetLineage(labels)
  
  # Check that the output is an empty list
  expect_type(lineage, "list")
  expect_length(lineage, 0)
})