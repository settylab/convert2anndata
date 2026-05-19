# Process Other Assays

This function processes assays other than the main assay in the
SingleCellExperiment object.

## Usage

``` r
process_other_assays(sce, assayName = NULL)
```

## Arguments

- sce:

  A SingleCellExperiment object.

- assayName:

  (Optional) The name of the assay used as the main data matrix. If not
  provided, all assays will be processed.

## Value

A list of processed assays.
