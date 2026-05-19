# Extract \`adata.raw\` as a CsparseMatrix

\`adata.raw\` (a snapshot of an AnnData often used by scanpy to retain
raw counts after normalisation) is a separate AnnData-like object with
its own \`.X\`, \`.var\`, and \`.var_names\`. It may have a different
number of features from the main AnnData object.

## Usage

``` r
extract_anndata_raw(adata, obs_names = NULL)
```

## Arguments

- adata:

  An AnnData object.

- obs_names:

  Character vector of cell names. Optional.

## Value

A CsparseMatrix of shape \`n_obs x n_raw_vars\`, or \`NULL\`.

## Details

This helper returns the raw matrix in AnnData orientation (cells x
features), with cell rownames from \`obs_names\` (the active AnnData)
and feature column names from \`adata\$raw\$var_names\`. Returns
\`NULL\` if \`adata\$raw\` is unset.
