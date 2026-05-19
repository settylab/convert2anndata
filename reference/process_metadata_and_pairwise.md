# Process Metadata and Pairwise Matrices

This function processes the metadata and pairwise matrices in a
\`SingleCellExperiment\` object. It organizes metadata into distinct
categories (\`uns\`, \`obsp\`, \`obsm\`, \`varp\`, \`varm\`) based on
their dimensions and relationship to the main data matrix (\`X\`).
Metadata matrices matching observation or variable dimensions are moved
into the corresponding slots, while remaining metadata is retained in
\`uns\`. Pairwise matrices for observations and variables are identified
and organized into \`obsp\` and \`varp\`, respectively.

## Usage

``` r
process_metadata_and_pairwise(sce, alt_exps_data, X, obsm = list())
```

## Arguments

- sce:

  A \`SingleCellExperiment\` object containing the data to be processed.

- alt_exps_data:

  A list of alternative experiment data to be added to the metadata.

- X:

  The main data matrix (\`anndata.X\`) for the AnnData object. This is
  used to validate dimensions of metadata and pairwise matrices.

- obsm:

  A list of observation matrices, which can be extended during
  processing.

## Value

A list containing:

- uns:

  A list of remaining metadata that was not categorized into other
  slots.

- obsp:

  A list of pairwise observation matrices.

- varp:

  A list of pairwise variable matrices.

- obsm:

  A list of observation matrices.

- varm:

  A list of variable matrices.
