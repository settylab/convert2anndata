library(testthat)
library(convert2anndata)
library(SingleCellExperiment)
library(Matrix)

test_that("process_main_assay works correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  X <- process_main_assay(sce, "counts")
  expect_equal(dim(X), c(2, 2))
})
