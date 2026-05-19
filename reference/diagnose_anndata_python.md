# Print a diagnostic block for the currently-active Python env

Convenience wrapper for \`check_anndata_python(verbose = TRUE, action =
"silent")\`. Always prints the diagnostic block; never raises. Returns
\`TRUE\` if the env is usable, \`FALSE\` otherwise.

## Usage

``` r
diagnose_anndata_python()
```

## Value

Logical, invisibly.
