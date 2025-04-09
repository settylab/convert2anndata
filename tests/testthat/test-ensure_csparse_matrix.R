test_that("ensure_csparse_matrix handles various matrix types correctly", {
  # Test with a dense matrix
  dense_mat <- matrix(1:9, 3, 3)
  result_dense <- ensure_csparse_matrix(dense_mat)
  expect_identical(dense_mat, result_dense)
  expect_false(methods::is(result_dense, "sparseMatrix"))
  
  # Test with dgCMatrix (already in CsparseMatrix format)
  dgc_mat <- Matrix::Matrix(1:9, 3, 3, sparse = TRUE)
  expect_true(methods::is(dgc_mat, "CsparseMatrix"))
  result_dgc <- ensure_csparse_matrix(dgc_mat)
  expect_identical(dgc_mat, result_dgc)
  expect_true(methods::is(result_dgc, "CsparseMatrix"))
  
  # Test with dgTMatrix (triplet format, not CsparseMatrix)
  ij <- expand.grid(i = 1:5, j = 1:5)
  ij <- ij[sample(nrow(ij), 10), ]
  dgt_mat <- Matrix::sparseMatrix(i = ij$i, j = ij$j, x = 1:10, dims = c(5, 5))
  expect_false(methods::is(dgt_mat, "CsparseMatrix"))
  result_dgt <- ensure_csparse_matrix(dgt_mat)
  expect_true(methods::is(result_dgt, "CsparseMatrix"))
  expect_equal(as.matrix(dgt_mat), as.matrix(result_dgt))
  
  # Test with a dgeMatrix (dense matrix from Matrix package)
  dge_mat <- Matrix::Matrix(matrix(1:9, 3, 3), sparse = FALSE)
  expect_false(methods::is(dge_mat, "sparseMatrix"))
  result_dge <- ensure_csparse_matrix(dge_mat)
  expect_identical(dge_mat, result_dge)
  
  # Test with a ddiMatrix (diagonal matrix)
  ddi_mat <- Matrix::Diagonal(3, 1:3)
  expect_false(methods::is(ddi_mat, "CsparseMatrix"))
  result_ddi <- ensure_csparse_matrix(ddi_mat)
  expect_true(methods::is(result_ddi, "CsparseMatrix"))
  expect_equal(as.matrix(ddi_mat), as.matrix(result_ddi))
})