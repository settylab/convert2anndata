# convert2anndata

[![Codecov test coverage](https://codecov.io/gh/settylab/convert2anndata/branch/main/graph/badge.svg)](https://codecov.io/gh/settylab/convert2anndata)

`convert2anndata` is an R package that provides functions to convert `SingleCellExperiment` and `Seurat` objects to `AnnData` format. It handles assays, dimensional reductions, metadata, pairing data such as cell-to-cell distances, and alternative experiments, ensuring comprehensive conversion.

## Installation

To install the `convert2anndata` package from GitHub through ssh, you can use the `remotes` package in R. If you don't have `remotes` installed, you can install it first:

```r
install.packages("remotes")
```

Then, install the necessary Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"))
```

Finally, install `convert2anndata` from GitHub:

```r
remotes::install_github("settylab/convert2anndata")
```

Alternatively, install from GitHub authenticating through ssh:

```r
remotes::install_git("git@github.com:settylab/convert2anndata.git")
```

### Installation with `renv`

You can set up the package and its dependencies in a project-specific environment using `renv`. This approach ensures that all dependencies are installed in a consistent environment.

First, install `renv` if you don't already have it:

```r
install.packages("renv")
```

Then, you can use the following steps to set up the environment and install all necessary packages, including Bioconductor packages, in one go:

```r
# Initialize renv in your project directory
renv::init(bare=TRUE)

# Install the necessary packages, including Bioconductor packages
renv::install(c(
    "bioc::SingleCellExperiment",
    "bioc::SummarizedExperiment",
    "bioc::S4Vectors"
))
    
renv::install("git@github.com:settylab/convert2anndata.git")

# Snapshot the environment
renv::snapshot(type="all")
```

### Note for Seurat Objects

If you want to convert Seurat objects, you will also need to install the Seurat package. Follow the installation instructions on the [Seurat website](https://satijalab.org/seurat/articles/install.html).


## Usage

### Command Line Usage

You can use the `convert2anndata` package from the command line to convert `SingleCellExperiment` or `Seurat` objects stored in RDS files to `AnnData` format (H5AD files). Here is an example of how to use it:

```sh
Rscript -e "convert2anndata::cli_convert()" -i /path/to/input_file.rds -o /path/to/output_file.h5ad
```

#### Command Line Options

- `-i`, `--input`: Path to the input RDS file containing the `SingleCellExperiment` or `Seurat` object. This option is required.
- `-o`, `--output`: Path to the output H5AD file. If not specified, the output path is derived by replacing the `.rds` extension of the input path with `.h5ad`.
- `-a`, `--assay`: The assay to use as the main matrix (`anndata.X`). Defaults to 'counts'.
- `-d`, `--disable-recursive-altExp`: Disable recursive recovery of `altExperiments` and discard them instead.

### R Usage

You can also use the `convert2anndata` package directly in R. Below are examples of how to convert `SingleCellExperiment` or `Seurat` objects to `AnnData` format within an R session.

#### Example: Converting a Seurat Object

```r
library(convert2anndata)

# Load a Seurat object
seurat_obj <- readRDS("/path/to/input_file.rds")

# Convert to SingleCellExperiment if necessary
sce <- convert_seurat_to_sce(seurat_obj)

# Convert to AnnData
ad <- convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE)

# Save the AnnData object
write_h5ad(ad, "/path/to/output_file.h5ad")
```

#### Example: Converting a SingleCellExperiment Object

```r
library(convert2anndata)

# Load a SingleCellExperiment object
sce <- readRDS("/path/to/input_file.rds")

# Convert to AnnData
ad <- convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE)

# Save the AnnData object
write_h5ad(ad, "/path/to/output_file.h5ad")
```

## Function Reference

### `convert_seurat_to_sce`

This function determines the class of a loaded object and converts it to a `SingleCellExperiment` object if necessary. It handles `Seurat` objects, updating old `Seurat` v2 objects if detected, and converts them to `SingleCellExperiment`. If the input object is already a `SingleCellExperiment`, it is returned as is. For other object types, an attempt is made to convert them to `SingleCellExperiment`.

### `convert_to_anndata`

This function converts a `SingleCellExperiment` (SCE) object to an `AnnData` object. It processes assays, dimensional reductions, and metadata, including alternative experiments (`altExps`). The main assay is used as the primary data matrix in the `AnnData` object.

### `cli_convert`

This function serves as a command line interface for the `convert2anndata` package. It parses command line arguments and calls the appropriate functions to convert a `SingleCellExperiment` or `Seurat` object to an `AnnData` object.

## License

This project is licensed under the GPL-3 License - see the [LICENSE](LICENSE) file for details.
