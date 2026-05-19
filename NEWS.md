# convert2anndata 0.3.0

## Bidirectional conversion

`convert2anndata` is now a bidirectional converter. In addition to the
existing SingleCellExperiment / Seurat → AnnData direction, it can
convert AnnData (`.h5ad`) objects into Seurat or SingleCellExperiment
objects:

* `convert_anndata_to_seurat()` — AnnData → Seurat.
* `convert_anndata_to_sce()` — AnnData → SingleCellExperiment.

Both accept either an in-memory AnnData object or a path to an `.h5ad`
file. AnnData components are mapped as follows: `X` / `layers` →
assays, `obsm` → dimensional reductions, `obs` / `var` → cell / feature
metadata, and `obsp` → Seurat graphs (for the Seurat target).

## Command-line interface

* `cli_convert()` now dispatches on the input file extension:
  `.rds` → `.h5ad` and `.h5ad` → `.rds`.
* The `.h5ad` → `.rds` direction supports `-t seurat` (default) or
  `-t sce` to choose the target object type.

## Python / reticulate setup

New helpers make the reticulate ↔ Python `anndata` wiring explicit and
easier to debug:

* `setup_anndata_python()` — resolve and activate a Python environment
  that has `anndata` installed.
* `check_anndata_python()` — fail fast with an actionable message when
  Python, `anndata`, or `numpy` are not reachable.
* `diagnose_anndata_python()` — print a full diagnostic of the active
  Python environment.

The path entry points to the converters call `check_anndata_python()`
automatically.

## Seurat 5 support and `use_raw`

* Conversions are robust across Seurat 4 / Seurat 5 layer and slot
  naming.
* `use_raw = "auto"` (the default) uses `adata.raw` as the counts
  matrix only when no requested `counts_layer` is present; `"always"`
  and `"never"` force the behaviour, and a logical value is still
  accepted for backward compatibility.

## New exported helpers

The full `extract_anndata_*` family is now exported for building custom
conversions: `extract_anndata_X()`, `extract_anndata_layers()`,
`extract_anndata_obsm()`, `extract_anndata_obsp()`, and
`extract_anndata_raw()`. Also new: `attach_reductions_seurat()`.

## Tests, evaluation harness, and validation

* About ten new test files cover the AnnData → Seurat / SCE direction,
  edge cases, sparse-matrix handling, `raw` / `obsp`, and realistic
  scanpy-style datasets. Tests that need Python `anndata` skip
  gracefully when it is unavailable.
* A standalone evaluation harness (`eval/`) exercises full roundtrips,
  edge cases, a real pbmc3k dataset, and a quantitative
  roundtrip-fidelity check.

### Validation results (0.3.0)

Validated with R 4.4.1, Seurat 5.4.0, SingleCellExperiment 1.28.1,
reticulate 1.45.0, and Python `anndata` 0.12.0:

* `R CMD check`: 0 errors, 0 warnings, 3 (benign) notes.
* `testthat`: 387 passing, 0 failing, 1 skipped (a CLI subprocess test
  that needs the package installed on the library path).
* Roundtrip fidelity: `X`, `obs`, `var`, `obsm`, and `obsp` roundtrip
  within numeric tolerance through the SCE-mediated path; `obsp` also
  roundtrips through the Seurat-mediated path.

### Known limitations

* The Seurat-mediated path (`AnnData → Seurat → ...`) cannot preserve
  arbitrary extra `layers` or `var` (feature) metadata, because a
  Seurat assay only has `counts` / `data` / `scale.data` slots and no
  feature-metadata store. Use the SCE-mediated path when verbatim
  layers or `var` columns matter.
* `convert_anndata_to_sce()` does not yet attach `obsp` matrices as
  `colPairs` (the Seurat target does attach them as graphs). This is a
  planned follow-up.

## Documentation

* `README` rewritten with a direction reference table, Python
  environment notes, and a troubleshooting section.
* `_pkgdown.yml` updated for bidirectional conversion with a grouped
  reference index.
* `DESCRIPTION` retitled and bumped to 0.3.0.
