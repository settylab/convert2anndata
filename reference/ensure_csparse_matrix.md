# Ensure Sparse Matrix is in a Format Compatible with AnnData

This function checks if a sparse matrix is in a format compatible with
AnnData (CsparseMatrix format, e.g., dgCMatrix) and converts it if
necessary.

## Usage

``` r
ensure_csparse_matrix(mat)
```

## Arguments

- mat:

  A matrix or Matrix object

## Value

The original matrix if already in a compatible format, or a converted
matrix if conversion was needed
