# tests/testthat/test_extract_counts_matrix.R

library(testthat)
library(Seurat)

test_that("extract_counts_matrix works with Assay (Seurat v3) objects", {
  skip_if_not_installed("Seurat")

  # Create an Assay object
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  assay <- CreateAssayObject(counts = counts_matrix)

  # Extract counts matrix
  counts <- extract_counts_matrix(assay)

  # Check that counts match
  expect_equal(
    as.matrix(counts),
    counts_matrix,
    check.attributes = FALSE
  )
})

test_that("extract_counts_matrix works with Assay5 (Seurat v5) objects", {
  skip_if_not_installed("Seurat")

  # Create an Assay5 object
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)
  assay <- seurat_obj@assays$RNA

  # Extract counts matrix
  counts <- extract_counts_matrix(assay)

  # Convert counts to matrix for comparison
  counts_matrix_extracted <- as.matrix(counts)

  # Check that counts match
  expect_equal(
    as.matrix(counts),
    counts_matrix,
    check.attributes = FALSE
  )
})
