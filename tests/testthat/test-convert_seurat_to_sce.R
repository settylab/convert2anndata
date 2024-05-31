library(testthat)
library(convert2anndata)
library(SingleCellExperiment)

test_that("convert_seurat_to_sce works with Seurat objects", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  seurat_obj <- Seurat::CreateSeuratObject(counts = matrix(1:4, ncol = 2))
  sce <- convert_seurat_to_sce(seurat_obj)
  expect_is(sce, "SingleCellExperiment")
})

test_that("convert_seurat_to_sce works with SingleCellExperiment objects", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  result <- convert_seurat_to_sce(sce)
  expect_is(result, "SingleCellExperiment")
})
