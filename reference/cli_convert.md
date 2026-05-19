# Command Line Interface for convert2anndata

Bidirectional CLI: dispatches on the input file extension. An \`.rds\`
input is treated as a Seurat or SingleCellExperiment object and
converted to an AnnData \`.h5ad\` file. An \`.h5ad\` input is read with
\`anndata::read_h5ad()\` and converted to a Seurat (default) or
SingleCellExperiment object saved to \`.rds\`.

## Usage

``` r
cli_convert()
```
