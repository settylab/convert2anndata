# Extract obsm embeddings from an AnnData object

Returns a named list of numeric matrices, one per \`adata\$obsm\` key.
Each matrix has rows in \`obs_names\` order. Non-numeric or
dimension-mismatched entries are skipped with a warning.

## Usage

``` r
extract_anndata_obsm(adata, obs_names = NULL)
```

## Arguments

- adata:

  An AnnData object.

- obs_names:

  Character vector of cell names; used both for sanity-checking row
  count and for setting rownames.

## Value

A named list of numeric matrices (cells x dim).
