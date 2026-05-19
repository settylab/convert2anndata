# Extract dict-like keys from an AnnData mapping (layers, obsm, uns, ...)

AnnData's mapping attributes (\`layers\`, \`obsm\`, \`uns\`, ...) expose
a \`.keys()\` method. Different combinations of \`reticulate\` and
\`anndata\` return that result in different shapes:

## Usage

``` r
anndata_mapping_keys(mapping)
```

## Arguments

- mapping:

  The mapping attribute (e.g. \`adata\$layers\`).

## Value

A character vector of keys, or \`character(0)\` on failure.

## Details

\- newer reticulate auto-converts the Python \`KeysView\` to an R
character vector (so \`as.character(builtins\$list(x\$keys()))\` runs
\`list()\` on a char-vec and splits each name into individual
characters); - older versions hand back a Python object that needs
\`builtins\$list()\` to materialise.

This helper handles both cases, plus the edge case where the mapping has
no \`.keys()\` method but \`names()\` works.
