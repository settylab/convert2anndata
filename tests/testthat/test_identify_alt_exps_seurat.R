# tests/testthat/test_identify_alt_exps_seurat.R

library(testthat)
library(Seurat)

test_that("identify_alt_exps_seurat identifies assays to move", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object with an additional assay
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  data <- CreateSeuratObject(counts = counts_matrix)

  # Add another assay
  counts_matrix2 <- matrix(21:40, nrow = 5, ncol = 4)
  rownames(counts_matrix2) <- paste0("Gene", 1:5)
  colnames(counts_matrix2) <- paste0("Cell", 1:4)
  assay2 <- CreateAssayObject(counts = counts_matrix2)
  data[["Assay2"]] <- assay2

  # Identify assays to move
  assays_to_move <- identify_alt_exps_seurat(data)

  # Should identify 'Assay2' to move
  expect_equal(assays_to_move, "Assay2")
})
