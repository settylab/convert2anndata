library(testthat)
library(convert2anndata)
library(SingleCellExperiment)
library(Matrix)

test_that("process_other_assays works correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2), logcounts = matrix(5:8, ncol = 2)))
  layers <- process_other_assays(sce, "counts")
  expect_equal(length(layers), 1)
  expect_equal(dim(layers$logcounts), c(2, 2))
})
