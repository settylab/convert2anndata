library(testthat)
library(convert2anndata)
library(SingleCellExperiment)

test_that("process_dimensional_reductions works correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  reducedDims(sce) <- SimpleList(PCA = matrix(1:4, ncol = 2))
  obsm <- process_dimensional_reductions(sce)
  expect_true("X_pca" %in% names(obsm))
  expect_equal(dim(obsm$X_pca), c(2, 2))
})
