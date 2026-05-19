# Extract obsp (cell x cell) matrices from an AnnData object

AnnData stores cell-pair matrices (e.g. nearest-neighbor graphs from
\`sc.pp.neighbors\`: \`connectivities\`, \`distances\`) in
\`adata.obsp\`. This helper returns them as a named list of sparse cell
x cell CsparseMatrix objects, ready to be attached as Seurat \`Graphs\`
or SCE \`colPairs\`.

## Usage

``` r
extract_anndata_obsp(adata, obs_names = NULL)
```

## Arguments

- adata:

  An AnnData object.

- obs_names:

  Character vector of cell names; used both for sanity- checking shape
  and for setting dimnames.

## Value

A named list of sparse matrices (\`n_obs x n_obs\`). Empty list if
\`adata\$obsp\` is empty / missing.
