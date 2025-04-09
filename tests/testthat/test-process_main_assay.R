test_that("process_main_assay converts non-CsparseMatrix to CsparseMatrix", {
  library(SingleCellExperiment)
  library(Matrix)
  
  # Create an assay with a triplet sparse matrix (dgTMatrix)
  triplet_matrix <- Matrix::sparseMatrix(
    i = c(1, 3, 5, 7),
    j = c(2, 4, 6, 8),
    x = 1:4,
    dims = c(10, 10)
  )
  expect_false(methods::is(triplet_matrix, "CsparseMatrix"))
  
  rownames(triplet_matrix) <- paste0("gene", 1:10)
  colnames(triplet_matrix) <- paste0("cell", 1:10)
  
  # Create SingleCellExperiment with triplet matrix
  sce <- SingleCellExperiment(list(counts = triplet_matrix))
  
  # Process the main assay
  X <- process_main_assay(sce, "counts")
  
  # Test that the result is a CsparseMatrix
  expect_true(methods::is(X, "CsparseMatrix"))
  
  # Test that data is preserved
  expected_values <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                       0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 3, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
  expect_equal(as.vector(as.matrix(X))[1:40], expected_values)
  
  # Test with a diagonal matrix
  diag_matrix <- Matrix::Diagonal(10, 1:10)
  expect_false(methods::is(diag_matrix, "CsparseMatrix"))
  
  rownames(diag_matrix) <- paste0("gene", 1:10)
  colnames(diag_matrix) <- paste0("cell", 1:10)
  
  # Create SingleCellExperiment with diagonal matrix
  sce2 <- SingleCellExperiment(list(counts = diag_matrix))
  
  # Process the main assay
  X2 <- process_main_assay(sce2, "counts")
  
  # Test that the result is a CsparseMatrix
  expect_true(methods::is(X2, "CsparseMatrix"))
  
  # Test that data is preserved
  expect_equal(diag(as.matrix(X2)), 1:10)
})