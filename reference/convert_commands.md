# Convert SeuratCommand objects into a serializable format

This function takes a list of SeuratCommand objects and converts them
into JSON-serializable lists, preserving command name, parameters,
timestamps, and the assay used.

## Usage

``` r
convert_commands(commands_list)
```

## Arguments

- commands_list:

  A list of SeuratCommand objects from Seurat metadata.

## Value

A list of simplified, serializable command information.
