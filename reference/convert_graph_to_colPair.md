# Convert Seurat Graph to colPair-Compatible Format

This function converts a graph (adjacency matrix or edge list) into a
format compatible with the \`colPair\` slot in a
\`SingleCellExperiment\` object. The graph must match the cell indices
or names in the \`SingleCellExperiment\` object.

## Usage

``` r
convert_graph_to_colPair(graph, colnames_sce)
```

## Arguments

- graph:

  A graph represented as either a square adjacency matrix, distance
  matrix, or an edge list (data frame or two-column matrix).

- colnames_sce:

  A character vector of cell names from the \`colData\` of the
  \`SingleCellExperiment\` object. Used to map the graph indices to the
  correct cell names.

## Value

A \`SelfHits\` object compatible with the \`colPair\` slot in
\`SingleCellExperiment\`.
