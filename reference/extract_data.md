# Extract Data from SingleCellExperiment

This function extracts data from the SingleCellExperiment object with
careful error handling.

## Usage

``` r
extract_data(data_fun, data_type, sce)
```

## Arguments

- data_fun:

  A function to extract data (e.g., colData, rowData).

- data_type:

  The type of data being extracted (e.g., "obs/colData", "var/rowData").

- sce:

  A SingleCellExperiment object.

## Value

A list containing the extracted data and internal columns.
