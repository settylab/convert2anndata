# Create Combined Assay in Seurat Object

Creates a new assay in the Seurat object using the combined counts
matrix.

## Usage

``` r
create_combined_assay(data, combined_counts, condition_labels)
```

## Arguments

- data:

  The Seurat object.

- combined_counts:

  The combined counts matrix.

- condition_labels:

  Labels indicating the origin of each cell.

## Value

The updated Seurat object.
