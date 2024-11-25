# tests/testthat/test_align_metadata_seurat.R

library(testthat)
library(Seurat)

test_that("align_metadata_seurat aligns metadata correctly", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object with meta.data
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  data <- CreateSeuratObject(counts = counts_matrix)
  data$cell_type <- c("A", "A", "B", "B")

  # Create combined_counts (reorder columns)
  combined_counts <- counts_matrix[, c(4, 3, 2, 1)]
  colnames(combined_counts) <- colnames(data)[c(4, 3, 2, 1)]

  # Align metadata
  metadata <- align_metadata_seurat(data, combined_counts)

  # Check that metadata is aligned
  expect_equal(
    unname(metadata$cell_type),
    unname(data$cell_type[c(4, 3, 2, 1)])
  )
  expect_equal(
    metadata$cell_type,
    data$cell_type[c(4, 3, 2, 1)],
    check.attributes = FALSE
  )
})

test_that("align_metadata_seurat throws error on mismatch", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object with meta.data
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  data <- CreateSeuratObject(counts = counts_matrix)
  data$cell_type <- c("A", "A", "B", "B")

  # Create combined_counts with non-matching cell names
  combined_counts <- counts_matrix[, 1:2]
  colnames(combined_counts) <- c("Cell5", "Cell6") # Cells not in data

  # Expect an error
  expect_error(
    align_metadata_seurat(data, combined_counts),
    "Mismatch between cell names in counts matrix and colData."
  )
})
