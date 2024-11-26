# tests/testthat/test_create_combined_assay.R

library(testthat)
library(Seurat)

test_that("create_combined_assay updates Seurat object correctly", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  data <- CreateSeuratObject(counts = counts_matrix)

  # Create combined_counts and condition_labels
  combined_counts <- counts_matrix
  rownames(combined_counts) <- rownames(counts_matrix)
  colnames(combined_counts) <- colnames(counts_matrix)
  condition_labels <- rep("TestCondition", ncol(combined_counts))

  # Call function
  data_updated <- create_combined_assay(data, combined_counts, condition_labels)

  # Check that new assay 'RNA_combined' is created
  expect_true("RNA_combined" %in% names(data_updated@assays))

  # Check that DefaultAssay is set to 'RNA_combined'
  expect_equal(DefaultAssay(data_updated), "RNA_combined")

  # Check that condition labels are added
  expect_true("condition" %in% colnames(data_updated@meta.data))
  expect_equal(
    unname(data_updated$condition),
    unname(factor(condition_labels))
  )
})
