# Extract obs/var names from an AnnData object

Cell and feature names on AnnData objects are Python pandas Indexes; R's
\`rownames()\`/\`colnames()\` do not return them reliably. This helper
coerces them to character vectors via Python's \`list()\` builtin,
falling back to \`names()\` if that path fails.

## Usage

``` r
anndata_names(adata)
```

## Arguments

- adata:

  An AnnData object (e.g. as returned by \`anndata::read_h5ad\`).

## Value

A list with elements \`obs_names\` (cells) and \`var_names\` (features),
each a character vector or \`NULL\` if extraction failed.
