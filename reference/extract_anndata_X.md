# Extract the primary count matrix from an AnnData object

Selects a layer from \`adata\$layers\` if \`layer\` is supplied and
present; otherwise falls back to \`adata\$X\`. The returned matrix is in
AnnData orientation (cells x features) and converted to CsparseMatrix.
Row and column names are set from \`obs_names\` / \`var_names\` when
provided.

## Usage

``` r
extract_anndata_X(adata, layer = NULL, obs_names = NULL, var_names = NULL)
```

## Arguments

- adata:

  An AnnData object.

- layer:

  Optional name of a layer in \`adata\$layers\` to use as the primary
  matrix. If \`NULL\` or missing, \`adata\$X\` is used.

- obs_names:

  Character vector of cell names (rows). Optional.

- var_names:

  Character vector of feature names (columns). Optional.

## Value

A CsparseMatrix of shape \`n_obs x n_vars\` (cells x features).
