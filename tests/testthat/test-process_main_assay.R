test_that("process_main_assay converts non-CsparseMatrix to CsparseMatrix", {
  library(SingleCellExperiment)
  library(Matrix)
  
  # Create an assay with a triplet sparse matrix (dgTMatrix). Matrix >= 1.5
  # defaults sparseMatrix() to dgCMatrix; force triplet via repr = "T".
  triplet_matrix <- Matrix::sparseMatrix(
    i = c(1, 3, 5, 7),
    j = c(2, 4, 6, 8),
    x = 1:4,
    dims = c(10, 10),
    repr = "T"
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
  
  # Test that data is preserved. process_main_assay transposes
  # (genes x cells -> cells x genes) before returning, so non-zeros at
  # M[1,2]=1, M[3,4]=2, M[5,6]=3, M[7,8]=4 become X[2,1]=1, X[4,3]=2,
  # X[6,5]=3, X[8,7]=4. Column-major linear indices: 2, 24, 46, 68 — only
  # the first two land in [1:40].
  expected_values <- c(
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,                # col 1 -> X[2,1]=1
    rep(0, 10),                                  # col 2
    0, 0, 0, 2, 0, 0, 0, 0, 0, 0,                # col 3 -> X[4,3]=2
    rep(0, 10)                                   # col 4
  )
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