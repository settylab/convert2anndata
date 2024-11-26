# tests/testthat/test_attach_alt_experiments_sce.R

library(testthat)
library(Seurat)
library(SingleCellExperiment)

test_that("attach_alt_experiments_sce attaches altExps correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")

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
  data[["Assay2_cell_count_mismatch"]] <- assay2

  # Create sce
  sce <- SingleCellExperiment(
    assays = list(counts = counts_matrix),
    colData = data@meta.data
  )

  # Attach altExps
  altExp_names <- "Assay2_cell_count_mismatch"
  sce <- attach_alt_experiments_sce(data, sce, altExp_names)

  # Check that altExp is attached
  expect_true("Assay2_cell_count_mismatch" %in% altExpNames(sce))

  # Check that altExp has correct counts
  alt_sce <- altExp(sce, "Assay2_cell_count_mismatch")
  expect_equal(
    as.matrix(assay(alt_sce, "counts")),
    counts_matrix2,
    check.attributes = FALSE
  )
})
