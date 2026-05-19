# Extract Pairwise Data from SingleCellExperiment

This function extracts pairwise data (colPairs or rowPairs) from a
SingleCellExperiment object and returns it as a list of sparse matrices.

## Usage

``` r
extract_pairs(pairs_function, pair_type, sce)
```

## Arguments

- pairs_function:

  A function to extract pairwise data (e.g., colPairs or rowPairs).

- pair_type:

  A string indicating whether the data is colPairing or rowPairing.

- sce:

  A SingleCellExperiment object from which to extract the pairwise data.

## Value

A list of pairwise data as sparse matrices.
