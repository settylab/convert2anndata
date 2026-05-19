# Default mapping from AnnData obsm keys to Seurat reduction (name, key) pairs

Provides the canonical Seurat naming for the obsm keys conventionally
produced by scanpy. Users who want extra mappings should pass a list
with the same shape into \`attach_reductions_seurat()\` (or its
callers).

## Usage

``` r
default_reduction_map()
```

## Value

A named list. Each name is an obsm key (e.g. \`X_pca\`); each element is
a list with \`name\` (Seurat reduction name, e.g. \`"pca"\`) and \`key\`
(column-name prefix, e.g. \`"PC\_"\`).
