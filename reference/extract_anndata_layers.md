# Extract layers from an AnnData object as a named list of matrices

Each layer is converted to CsparseMatrix and given dimnames from
\`obs_names\` / \`var_names\`. The orientation is preserved as
AnnData-native (cells x features); callers that need genes x cells (e.g.
Seurat) should transpose.

## Usage

``` r
extract_anndata_layers(
  adata,
  obs_names = NULL,
  var_names = NULL,
  exclude = character(0)
)
```

## Arguments

- adata:

  An AnnData object.

- obs_names:

  Character vector of cell names.

- var_names:

  Character vector of feature names.

- exclude:

  Character vector of layer names to skip (e.g. one already used as the
  primary \`X\`).

## Value

A named list of CsparseMatrix objects (cells x features). Empty list if
the AnnData has no layers.
