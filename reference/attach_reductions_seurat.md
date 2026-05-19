# Attach obsm embeddings to a Seurat object as dimensional reductions

Uses an obsm-key -\> (reduction-name, key-prefix) mapping. Keys not in
the mapping fall back to a derived name: leading \`X\_\` stripped,
lowercased, with a sanitized key prefix.

## Usage

``` r
attach_reductions_seurat(
  seurat_obj,
  obsm,
  assay = "RNA",
  reduction_map = list()
)
```

## Arguments

- seurat_obj:

  A Seurat object to attach reductions to.

- obsm:

  Named list of numeric matrices (cells x dim), as returned by
  \`extract_anndata_obsm()\`.

- assay:

  Assay name to associate the reductions with. Defaults to "RNA".

- reduction_map:

  Named list overriding or extending \`default_reduction_map()\`. Each
  entry maps an obsm key (e.g. \`"X_pca"\`) to a list with \`name\` and
  \`key\`. The supplied map is merged on top of the defaults; pass an
  empty list to keep defaults, or pass an entry with \`name = NA\` to
  disable a default.

## Value

The Seurat object with reductions attached.
