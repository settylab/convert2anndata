# convert2anndata

[![Codecov test
coverage](https://codecov.io/gh/settylab/convert2anndata/branch/main/graph/badge.svg)](https://codecov.io/gh/settylab/convert2anndata)

`convert2anndata` is an R package for **bidirectional conversion**
between `AnnData` (the canonical Python single-cell format) and either
`SingleCellExperiment` or `Seurat` objects. It handles split layers
(Seurat), assays, dimensional reductions (`obsm` ↔︎ `reducedDims` /
Seurat reductions), metadata (`obs`/`var` ↔︎ `colData`/`rowData`),
layers, and alternative experiments, aiming for a faithful roundtrip. If
you encounter any issues or notice incomplete conversions, please feel
free to report them on our [GitHub issue
tracker](https://github.com/settylab/convert2anndata/issues) to help us
continuously improve.

### Direction reference

| From → To              | Function                                                                                                                                                     |
|------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Seurat / SCE → AnnData | `convert_to_anndata(sce, ...)` (Seurat first via [`convert_seurat_to_sce()`](https://settylab.github.io/convert2anndata/reference/convert_seurat_to_sce.md)) |
| AnnData → Seurat       | `convert_anndata_to_seurat(adata, ...)`                                                                                                                      |
| AnnData → SCE          | `convert_anndata_to_sce(adata, ...)`                                                                                                                         |

## Installation

To install the `convert2anndata` package from GitHub through ssh, you
can use the `remotes` package in R. If you don’t have `remotes`
installed, you can install it first:

``` r
install.packages("remotes")
```

Then, install the necessary Bioconductor packages:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"))
```

Finally, install `convert2anndata` from GitHub:

``` r
remotes::install_github("settylab/convert2anndata")
```

Alternatively, install from GitHub authenticating through ssh:

``` r
remotes::install_git("git@github.com:settylab/convert2anndata.git")
```

### Installation with `renv`

You can set up the package and its dependencies in a project-specific
environment using `renv`. This approach ensures that all dependencies
are installed in a consistent environment.

First, install `renv` if you don’t already have it:

``` r
install.packages("renv")
```

Then, you can use the following steps to set up the environment and
install all necessary packages, including Bioconductor packages, in one
go:

``` r
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

If you want to convert Seurat objects, you will also need to install the
Seurat package. Follow the installation instructions on the [Seurat
website](https://satijalab.org/seurat/articles/install.html).

### Make an alias

Consider making an alias for the command line tool, e.g., with

``` bash
alias c2a='Rscript -e "convert2anndata::cli_convert()"'
echo 'alias c2a="Rscript -e \"convert2anndata::cli_convert()\""' >> ~/.bashrc
```

Now you can use the command line toole, explained under
`Command Line Usage` below, just by typing, e.g., `c2a -h`.

## Usage

### Command Line Usage

The CLI dispatches on the input file extension. `.rds` input is treated
as a Seurat or SingleCellExperiment object and converted to `.h5ad`;
`.h5ad` input is read with
[`anndata::read_h5ad()`](https://anndata.dynverse.org/reference/read_h5ad.html)
and converted to a Seurat (default) or SCE object saved to `.rds`.

``` sh
# RDS -> H5AD
Rscript -e "convert2anndata::cli_convert()" -i /path/to/input.rds -o /path/to/output.h5ad

# H5AD -> RDS (Seurat by default)
Rscript -e "convert2anndata::cli_convert()" -i /path/to/input.h5ad -o /path/to/output.rds

# H5AD -> RDS (SingleCellExperiment)
Rscript -e "convert2anndata::cli_convert()" -i /path/to/input.h5ad -t sce
```

If you set up an alias, as suggested in the `Installation` section, then
you can also conviniently run

``` sh
c2a -i /path/to/input_file.rds -o /path/to/output_file.h5ad
```

#### Command Line Options

- `-i`, `--input`: Path to the input file (`.rds` or `.h5ad`). Required.
- `-o`, `--output`: Path to the output file. If not specified, the
  output path is derived by swapping the extension (`.rds` ↔︎ `.h5ad`).
- `-a`, `--assay`: For `.rds → .h5ad`, the assay to use as `anndata.X`.
  For `.h5ad → .rds`, the layer name to use as the counts assay.
  Defaults to `counts`.
- `-d`, `--disable-recursive-altExp`: Disable recursive recovery of
  `altExperiments` and discard them instead. Applies to `.rds → .h5ad`.
- `-t`, `--target`: For `.h5ad → .rds` only: target object type,
  `seurat` (default) or `sce`.
- `-h`, `--help`: Show a help message and exit.

### R Usage

You can also use the `convert2anndata` package directly in R. Below are
examples of how to convert `SingleCellExperiment` or `Seurat` objects to
`AnnData` format within an R session.

#### Seurat / SCE → AnnData

``` r
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

#### AnnData → Seurat

``` r
library(convert2anndata)
library(anndata)

# Read an .h5ad file as an AnnData R6 object
adata <- read_h5ad("/path/to/input.h5ad")

# Convert to Seurat. counts_layer selects which layer becomes the counts
# assay; if the layer is missing, adata$X is used.
seurat_obj <- convert_anndata_to_seurat(adata, counts_layer = "counts")

saveRDS(seurat_obj, "/path/to/output.rds")
```

#### AnnData → SingleCellExperiment

``` r
adata <- read_h5ad("/path/to/input.h5ad")
sce <- convert_anndata_to_sce(adata, counts_layer = "counts")
```

Find the function documentation in the [reference
manual](https://settylab.github.io/convert2anndata/reference/) or
retrive the documentation through
[`?convert_to_anndata`](https://settylab.github.io/convert2anndata/reference/convert_to_anndata.md)
for any of functions.

## Python environment

`convert2anndata` is a thin R-side wrapper around the Python `anndata`
library, accessed through `reticulate`. **You must have a Python
interpreter that has `anndata` (and a compatible `numpy`) installed
before any conversion call will work.** The package does not bundle
Python.

The simplest setup uses the bundled installer from the `anndata` R
package, which provisions a managed Python venv for you:

``` r
anndata::install_anndata()
```

If you already have a Python environment with `anndata` installed
(e.g. a conda/micromamba env, or a `uv` venv), point reticulate at it.
The package looks at, in order:

1.  an explicit `conda_env` argument (path entrypoints only),
2.  `Sys.getenv("RETICULATE_PYTHON")`,
3.  `Sys.getenv("CONDA_PREFIX")`.

[`setup_anndata_python()`](https://settylab.github.io/convert2anndata/reference/setup_anndata_python.md)
is a small helper that runs the same resolution chain and emits a
warning instead of a hard error on misconfiguration:

``` r
convert2anndata::setup_anndata_python("my-anndata-env")
# or:
Sys.setenv(RETICULATE_PYTHON = "/path/to/python")
```

Use
[`check_anndata_python()`](https://settylab.github.io/convert2anndata/reference/check_anndata_python.md)
at the start of a script to fail fast with an actionable message if
Python, the `anndata` module, or numpy aren’t reachable. The path
entrypoints (`convert_anndata_to_seurat("file.h5ad")` and
`convert_anndata_to_sce("file.h5ad")`) call this automatically.

### Troubleshooting

**`No Python interpreter available to reticulate.`** Reticulate could
not find any Python. Set `RETICULATE_PYTHON` or run
[`anndata::install_anndata()`](https://anndata.dynverse.org/reference/install_anndata.html).

**`Python module 'anndata' is not importable.`** Reticulate found a
Python, but `anndata` isn’t installed in it. Either install into the
active env (`reticulate::py_install("anndata")`) or point
`RETICULATE_PYTHON` at a different one.

**`ImportError: Error importing numpy: you should not try to import numpy from its source directory`**
Almost always caused by a polluted `PYTHONPATH`. Common on HPC where
loading an R module also exports site-packages paths for a *different*
Python version (e.g. `/app/software/SciPy-bundle/.../python3.11/...`).
The fix is to clear the variable before launching R:

``` sh
unset PYTHONPATH
Rscript your-script.R
```

Or from inside R, **before any reticulate import**:

``` r
Sys.unsetenv("PYTHONPATH")
```

(after which a fresh R session is safest).

**`use_condaenv("missing-env")` silently “succeeds” then explodes
later.** With `required = FALSE` (the default), reticulate records the
preference without validating it. Either pass `required = TRUE` or rely
on
[`check_anndata_python()`](https://settylab.github.io/convert2anndata/reference/check_anndata_python.md)
to fail fast at the next call.

**Differences between conda envs:** `anndata` 0.8 / 0.9 / 0.12 differ
slightly. The package is tested against 0.8.x. If you hit a
method-not-found error from inside conversion, check
[`reticulate::py_config()`](https://rstudio.github.io/reticulate/reference/py_config.html)
and the Python `anndata.__version__`.

## License

This project is licensed under the GPL-3 License - see the
[LICENSE](https://settylab.github.io/convert2anndata/LICENSE) file for
details.
