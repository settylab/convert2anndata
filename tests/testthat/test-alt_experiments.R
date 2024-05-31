library(testthat)
library(convert2anndata)
library(SingleCellExperiment)

test_that("process_alt_experiments works correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  alt_sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(5:8, ncol = 2)))
  altExp(sce, "alt1") <- alt_sce
  result <- process_alt_experiments(sce, "counts", TRUE)
  expect_is(result$sce, "SingleCellExperiment")
  expect_true("alt1" %in% names(result$alt_exps))
})
