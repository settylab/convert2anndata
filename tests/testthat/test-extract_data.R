library(testthat)
library(convert2anndata)
library(SingleCellExperiment)

test_that("extract_data works correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  obs_data <- extract_data(colData, "obs/colData", sce)
  expect_is(obs_data$data, "data.frame")
})
