# Identify Assays to Move to altExps in Seurat Object

Identifies assays with mismatching cells or features and prepares them
to be moved to \`altExps\`.

## Usage

``` r
identify_alt_exps_seurat(data)
```

## Arguments

- data:

  The Seurat object.

## Value

A vector of assay names to be moved to \`altExps\`.
