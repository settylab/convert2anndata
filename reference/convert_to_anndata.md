# Convert SingleCellExperiment to AnnData

This function converts a SingleCellExperiment (SCE) object to an AnnData
object. It processes assays, dimensional reductions, and metadata,
including alternative experiments (altExps). The main assay is used as
the primary data matrix in the AnnData object.

## Usage

``` r
convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE)
```

## Arguments

- sce:

  A SingleCellExperiment object to be converted.

- assayName:

  The name of the assay to use as the main data matrix in AnnData.
  Defaults to "counts".

- useAltExp:

  Logical indicating whether to process and include alternative
  experiments (altExps). Defaults to TRUE.

## Value

An AnnData object containing the data from the SingleCellExperiment
object.

## Details

The function first prints a summary of the input SingleCellExperiment
object and checks for alternative experiments (altExps). If altExps are
found, they are processed recursively and stored in the \`uns\` slot of
the AnnData object, unless the \`useAltExp\` parameter is set to FALSE.
The specified assay is used as the primary data matrix (anndata.X), and
all other assays are stored in the \`layers\` slot. Dimensional
reductions are stored in the \`obsm\` slot, and metadata is organized
into \`obs\`, \`var\`, \`obsp\`, \`varp\`, and \`uns\` slots as
appropriate. The function includes detailed error handling to ensure
that all data is correctly transferred to the AnnData object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(SingleCellExperiment)
library(anndata)
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
ad <- convert_to_anndata(sce)
} # }
```
