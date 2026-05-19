# Check that reticulate has a working Python with the \`anndata\` module

Reticulate-based packages typically fail late and cryptically when the
Python side isn't wired up correctly: the user gets a numpy import error
from deep inside Python, far from the call site, with no hint as to what
to fix. This helper runs a layered check at the start of a conversion
call and stops with a single message that says exactly what is wrong and
how to fix it.

## Usage

``` r
check_anndata_python(action = c("stop", "warn", "silent"), verbose = FALSE)
```

## Arguments

- action:

  What to do on failure: \`"stop"\` (default) raises an error;
  \`"warn"\` emits a warning and returns \`FALSE\`; \`"silent"\` returns
  \`FALSE\` quietly.

- verbose:

  Logical. When \`TRUE\`, prints a full diagnostic block on success
  showing Python path, version, anndata version, numpy version,
  \`RETICULATE_PYTHON\` and \`CONDA_PREFIX\`. Useful as a one-stop "what
  does my env actually look like" call.

## Value

\`TRUE\` (invisibly) when all checks pass, \`FALSE\` when any check
fails and \`action != "stop"\`.

## Details

Checks performed, in order, with stop-at-first-failure semantics:

1\. \*\*A Python interpreter is reachable.\*\* If
\`reticulate::py_available()\` is \`FALSE\` even after initialisation,
no Python could be located. 2. \*\*\`PYTHONPATH\` does not point at a
different Python's site-packages.\*\* A common HPC failure mode is that
loading an R module exports paths like
\`/app/.../Python-3.11/site-packages\` while reticulate has selected a
3.8 / 3.12 interpreter, producing the cryptic "Error importing numpy:
you should not try to import numpy from its source directory" error. We
detect this \*before\* it fires. 3. \*\*The Python \`anndata\` module is
importable.\*\* 4. \*\*\`numpy\` imports without complaint\*\* (catches
any remaining ABI mismatches). 5. \*\*An end-to-end smoke probe\*\*:
round-trip a 2x2 AnnData object through reticulate. This is the only way
to be sure the env is actually usable for conversion work.

Each failure produces a tailored error message including the Python
interpreter that was selected and a one-line suggested fix.

## Examples

``` r
if (FALSE) { # \dontrun{
check_anndata_python()                    # stops with a helpful message
check_anndata_python(verbose = TRUE)      # also prints env diagnostics
if (check_anndata_python(action = "silent")) {
  adata <- anndata::read_h5ad(path)
}
} # }
```
