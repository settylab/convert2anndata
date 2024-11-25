# convert2anndata

[![Codecov test coverage](https://codecov.io/gh/settylab/convert2anndata/branch/main/graph/badge.svg)](https://codecov.io/gh/settylab/convert2anndata)

`convert2anndata` is an R package designed to seamlessly convert `SingleCellExperiment` and `Seurat` objects into the `AnnData` format, widely used in single-cell data analysis. The package supports the conversion of split layers (Seurat), assays, dimensional reductions, metadata, cell-to-cell pairing data (e.g., distances), and alternative experiments, ensuring a comprehensive transfer of information. If you encounter any issues or notice incomplete conversions, please feel free to report them on our [GitHub issue tracker](https://github.com/settylab/convert2anndata/issues) to help us continuously improve.


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

### Make an alias

Consider making an alias for the command line tool, e.g., with

```bash
alias c2a='Rscript -e "convert2anndata::cli_convert()"'
echo 'alias c2a="Rscript -e \"convert2anndata::cli_convert()\""' >> ~/.bashrc
```

Now you can use the command line toole, explained under `Command Line Usage` below, just by typing, e.g., `c2a -h`.

## Usage

### Command Line Usage

You can use the `convert2anndata` package from the command line to convert `SingleCellExperiment` or `Seurat` objects stored in RDS files to `AnnData` format (H5AD files). Here is an example of how to use it:

```sh
Rscript -e "convert2anndata::cli_convert()" -i /path/to/input_file.rds -o /path/to/output_file.h5ad
```

If you set up an alias, as suggested in the `Installation` section, then you can also conviniently run

```sh
c2a -i /path/to/input_file.rds -o /path/to/output_file.h5ad
```

#### Command Line Options

- `-i`, `--input`: Path to the input RDS file containing the `SingleCellExperiment` or `Seurat` object. This option is required.
- `-o`, `--output`: Path to the output H5AD file. If not specified, the output path is derived by replacing the `.rds` extension of the input path with `.h5ad`.
- `-a`, `--assay`: The assay to use as the main matrix (`anndata.X`). Defaults to 'counts'.
- `-d`, `--disable-recursive-altExp`: Disable recursive recovery of `altExperiments` and discard them instead.
- `-h`, `--help`: Show a help massage and exit.

### R Usage

You can also use the `convert2anndata` package directly in R. Below are examples of how to convert `SingleCellExperiment` or `Seurat` objects to `AnnData` format within an R session.

#### Example

```r
library(convert2anndata)
library(anndata)

# Load a Seurat object
seurat_obj <- readRDS("/path/to/input_file.rds")

# Convert to SingleCellExperiment if necessary
sce <- convert_seurat_to_sce(seurat_obj)

# Convert to AnnData
ad <- convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE)

# Save the AnnData object
write_h5ad(ad, "/path/to/output_file.h5ad")
```

Find the function documentation in the [reference manual](https://settylab.github.io/convert2anndata/reference/)
or retrive the documentation through `?convert_to_anndata` for any of functions.

## License

This project is licensed under the GPL-3 License - see the [LICENSE](LICENSE) file for details.
