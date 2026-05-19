# Process Alternative Experiments

This function processes alternative experiments (altExps) in the
SingleCellExperiment object. If altExps are found, they are processed
recursively and stored in the \`uns\` slot of the AnnData object.

## Usage

``` r
process_alt_experiments(sce, assayName, useAltExp)
```

## Arguments

- sce:

  A SingleCellExperiment object.

- assayName:

  The name of the assay to use as the main data matrix in AnnData.

- useAltExp:

  Logical indicating whether to process and include alternative
  experiments (altExps).

## Value

A list containing the updated SingleCellExperiment object and the
altExps data.
