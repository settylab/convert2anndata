library(testthat)
library(convert2anndata)
library(SingleCellExperiment)
library(anndata)
library(R6)

test_that("convert_to_anndata works correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  tryCatch({
    ad <- convert_to_anndata(sce)
    expect_true(R6::is.R6(ad)) # Ensure it is an R6 object
    expect_equal(nrow(ad$X), 2)
    expect_equal(ncol(ad$X), 2)
    if (!is.null(ad$obs)) {
      expect_equal(nrow(ad$obs), 2)
    }
    if (!is.null(ad$var)) {
      expect_equal(nrow(ad$var), 2)
    }
  }, error = function(e) {
    print(reticulate::py_last_error())
    stop(e)
  })
})

test_that("convert_to_anndata handles altExps correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  alt_sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(5:8, ncol = 2)))
  altExp(sce, "alt1") <- alt_sce
  tryCatch({
    ad <- convert_to_anndata(sce, useAltExp = TRUE)
    expect_true("altExperiments" %in% names(ad$uns))
    expect_true(R6::is.R6(ad$uns$altExperiments$alt1)) # Ensure it is an R6 object
    if (!is.null(ad$obs)) {
      expect_equal(nrow(ad$obs), 2)
    }
    if (!is.null(ad$var)) {
      expect_equal(nrow(ad$var), 2)
    }
  }, error = function(e) {
    print(reticulate::py_last_error())
    stop(e)
  })
})

test_that("convert_to_anndata works with Seurat objects", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  seurat_obj <- Seurat::CreateSeuratObject(counts = matrix(1:4, ncol = 2))
  sce <- convert_seurat_to_sce(seurat_obj)
  ad <- convert_to_anndata(sce)
  expect_true(R6::is.R6(ad)) # Ensure it is an R6 object
  expect_equal(nrow(ad$X), 2)
  expect_equal(ncol(ad$X), 2)
  if (!is.null(ad$obs)) {
    expect_equal(nrow(ad$obs), 2)
  }
  if (!is.null(ad$var)) {
    expect_equal(nrow(ad$var), 2)
  }
})
