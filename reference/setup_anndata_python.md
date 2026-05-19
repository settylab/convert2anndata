# Configure the Python interpreter used by reticulate for AnnData I/O

Call this once at the start of an R session to point reticulate at a
Python interpreter that has the \`anndata\` module installed. The
resolution order is:

## Usage

``` r
setup_anndata_python(conda_env = NULL, required = FALSE, validate = TRUE)
```

## Arguments

- conda_env:

  Optional. A conda environment name (resolved against common locations:
  micromamba, miniconda3, conda) or an absolute path to an environment
  prefix.

- required:

  Logical; passed to the underlying reticulate calls. Defaults to
  \`FALSE\` so that misconfiguration emits a warning rather than
  aborting.

- validate:

  Logical. After setting the interpreter preference, run
  \`check_anndata_python(action = "warn")\` to confirm the env actually
  imports \`anndata\`. Defaults to \`TRUE\`. With \`required = FALSE\`
  reticulate's own validation is a no-op, so this is the only protection
  against silent typos in \`conda_env\`.

## Value

Invisibly returns a single-element character vector describing which
path was used: one of \`"conda_env"\`, \`"conda_env_resolved"\`,
\`"reticulate_python"\`, \`"conda_prefix"\`, or \`"none"\`.

## Details

1\. \`conda_env\` (if supplied) — first tried as a conda environment
name, then as a path to an environment prefix, then as a path under
\`~/micromamba/envs/\<name\>\` or \`~/miniconda3/envs/\<name\>\`. 2.
\`RETICULATE_PYTHON\` environment variable — used directly as a Python
binary path if it exists. 3. \`CONDA_PREFIX\` environment variable —
treated as the active conda/micromamba environment prefix. 4. No-op
(reticulate falls back to its own discovery).

Failures emit a warning and continue rather than throwing, so users who
already have Python wired up don't have to special-case this call.

## Examples

``` r
if (FALSE) { # \dontrun{
setup_anndata_python("my-anndata-env")
adata <- anndata::read_h5ad("pbmc.h5ad")
} # }
```
