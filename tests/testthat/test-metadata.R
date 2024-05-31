library(testthat)
library(convert2anndata)
library(SingleCellExperiment)

test_that("process_metadata_and_pairwise works correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  X <- matrix(1:4, ncol = 2)
  uns_data <- process_metadata_and_pairwise(sce, list(), X)
  expect_is(uns_data$uns, "list")
})
