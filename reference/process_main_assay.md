# Process Main Assay

This function processes the main assay in the SingleCellExperiment
object to be used as the primary data matrix in the AnnData object.

## Usage

``` r
process_main_assay(sce, assayName)
```

## Arguments

- sce:

  A SingleCellExperiment object.

- assayName:

  The name of the assay to use as the main data matrix in AnnData.

## Value

The main data matrix (anndata.X) for the AnnData object.
