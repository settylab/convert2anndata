library(testthat)
library(Seurat)

test_that("update_seurat_object leaves Seurat v3+ object unchanged", {
  skip_if_not_installed("Seurat")

  # Create a Seurat v3 object
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  seurat_obj_v3 <- CreateSeuratObject(counts = counts_matrix)
  seurat_obj_v3@version <- package_version("3.1.0")

  # Update the Seurat object
  updated_obj <- update_seurat_object(seurat_obj_v3)

  # Check that the object is unchanged
  expect_identical(updated_obj, seurat_obj_v3)
})
