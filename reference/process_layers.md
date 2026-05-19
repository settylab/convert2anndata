# Process Layers in Seurat Assay

Combines split layers in a Seurat Assay5 object, stores them as
non-split layers, and removes the split layers.

## Usage

``` r
process_layers(assay_object)
```

## Arguments

- assay_object:

  The Seurat Assay5 object to process.

## Value

A list with split names per layer and updated assay object.
