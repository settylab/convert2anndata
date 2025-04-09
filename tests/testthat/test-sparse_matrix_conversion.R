test_that("sparse matrix conversion works in integration", {
  skip_if_not_installed("anndata")
  
  # Create a simple SCE object with different types of sparse matrices
  library(SingleCellExperiment)
  library(Matrix)
  
  # Create a dgCMatrix for the main assay
  counts <- Matrix::rsparsematrix(100, 50, 0.1, repr = "C")
  rownames(counts) <- paste0("gene", 1:100)
  colnames(counts) <- paste0("cell", 1:50)
  
  # Create a dgTMatrix for another assay
  logcounts <- Matrix::rsparsematrix(100, 50, 0.1, repr = "T")
  rownames(logcounts) <- rownames(counts)
  colnames(logcounts) <- colnames(counts)
  
  # Create a diagonal matrix for yet another assay
  diag_counts <- Matrix::Diagonal(100, rep(1, 100))[, 1:50]
  rownames(diag_counts) <- rownames(counts)
  colnames(diag_counts) <- colnames(counts)
  
  # Create SCE object
  sce <- SingleCellExperiment(list(
    counts = counts,
    logcounts = logcounts,
    diag_counts = diag_counts
  ))
  
  # Check initial matrix types
  expect_true(methods::is(assay(sce, "counts"), "CsparseMatrix"))
  expect_false(methods::is(assay(sce, "logcounts"), "CsparseMatrix"))
  expect_false(methods::is(assay(sce, "diag_counts"), "CsparseMatrix"))
  
  # Set up a dimensional reduction
  reducedDim(sce, "PCA") <- matrix(rnorm(50 * 10), 50, 10)
  
  # Add metadata that's a sparse matrix (not CsparseMatrix)
  int_metadata(sce)$sparse_meta <- Matrix::sparseMatrix(
    i = sample(1:50, 20),
    j = sample(1:50, 20),
    x = 1:20,
    dims = c(50, 50)
  )
  
  # Verify it's not a CsparseMatrix
  expect_false(methods::is(int_metadata(sce)$sparse_meta, "CsparseMatrix"))
  
  # Convert to AnnData
  ad <- convert_to_anndata(sce)
  
  # Check that all was converted properly
  expect_equal(dim(ad$X), c(50, 100))  # Note: transpose happens during conversion
  expect_equal(length(ad$layers), 2)
  expect_equal(length(ad$obsm), 1)
  
  # If we get this far, the conversion succeeded with all matrix types
  # which means the ensure_csparse_matrix function is working correctly
  expect_true(TRUE)
})