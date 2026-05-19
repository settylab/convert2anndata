# Preserve Dimensional Reductions from Seurat Object

Updates and transfers dimensional reductions from the Seurat object to
the SingleCellExperiment object by subsetting and reordering the cell
embeddings to match the cells in the SCE. For reductions that cover only
a subset of cells, missing rows are padded with NA.

## Usage

``` r
preserve_reductions_seurat(data, sce)
```

## Arguments

- data:

  The Seurat object.

- sce:

  The SingleCellExperiment object.

## Value

The SingleCellExperiment object with updated dimensional reductions.
