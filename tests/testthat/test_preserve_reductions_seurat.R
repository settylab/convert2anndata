library(testthat)
library(Seurat)
library(SingleCellExperiment)

test_that("preserve_reductions_seurat transfers reductions correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")

  # Create a Seurat object with reductions
  counts_matrix <- matrix(rpois(100, lambda = 10), nrow = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:10)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)

  # Add PCA reduction
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 5)

  # Create sce
  sce <- SingleCellExperiment(
    assays = list(counts = counts_matrix),
    colData = seurat_obj@meta.data
  )

  # Preserve reductions
  sce <- preserve_reductions_seurat(seurat_obj, sce)

  # Check that reductions are transferred
  expect_true("pca" %in% reducedDimNames(sce))
  expect_equal(
    reducedDim(sce, "pca"),
    Embeddings(seurat_obj, "pca")[colnames(seurat_obj), ]
  )
})
