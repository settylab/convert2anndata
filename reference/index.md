# Package index

## Bidirectional conversion

Main entry points. Converters accept either an in-memory object or a
path to an .h5ad / .rds file; cli_convert() is the command-line
dispatcher that picks a direction from the file extension.

- [`convert_anndata_to_seurat()`](https://settylab.github.io/convert2anndata/reference/convert_anndata_to_seurat.md)
  : Convert an AnnData object (or a .h5ad file) to a Seurat object
- [`convert_anndata_to_sce()`](https://settylab.github.io/convert2anndata/reference/convert_anndata_to_sce.md)
  : Convert an AnnData object (or a .h5ad file) to a
  SingleCellExperiment
- [`convert_to_anndata()`](https://settylab.github.io/convert2anndata/reference/convert_to_anndata.md)
  : Convert SingleCellExperiment to AnnData
- [`convert_seurat_to_sce()`](https://settylab.github.io/convert2anndata/reference/convert_seurat_to_sce.md)
  : Convert Seurat or Other Object to SingleCellExperiment
- [`cli_convert()`](https://settylab.github.io/convert2anndata/reference/cli_convert.md)
  : Command Line Interface for convert2anndata

## Python environment

Helpers for wiring reticulate to a Python interpreter that has the
anndata module installed, and for diagnosing misconfiguration.

- [`setup_anndata_python()`](https://settylab.github.io/convert2anndata/reference/setup_anndata_python.md)
  : Configure the Python interpreter used by reticulate for AnnData I/O
- [`check_anndata_python()`](https://settylab.github.io/convert2anndata/reference/check_anndata_python.md)
  : Check that reticulate has a working Python with the \`anndata\`
  module
- [`diagnose_anndata_python()`](https://settylab.github.io/convert2anndata/reference/diagnose_anndata_python.md)
  : Print a diagnostic block for the currently-active Python env

## AnnData component extractors

Lower-level accessors that pull individual pieces out of an AnnData
object as R-native matrices / data frames. Useful when building a custom
conversion.

- [`extract_anndata_X()`](https://settylab.github.io/convert2anndata/reference/extract_anndata_X.md)
  : Extract the primary count matrix from an AnnData object
- [`extract_anndata_layers()`](https://settylab.github.io/convert2anndata/reference/extract_anndata_layers.md)
  : Extract layers from an AnnData object as a named list of matrices
- [`extract_anndata_obsm()`](https://settylab.github.io/convert2anndata/reference/extract_anndata_obsm.md)
  : Extract obsm embeddings from an AnnData object
- [`extract_anndata_obsp()`](https://settylab.github.io/convert2anndata/reference/extract_anndata_obsp.md)
  : Extract obsp (cell x cell) matrices from an AnnData object
- [`extract_anndata_raw()`](https://settylab.github.io/convert2anndata/reference/extract_anndata_raw.md)
  : Extract \`adata.raw\` as a CsparseMatrix
- [`anndata_names()`](https://settylab.github.io/convert2anndata/reference/anndata_names.md)
  : Extract obs/var names from an AnnData object

## Attaching components to R objects

- [`attach_reductions_seurat()`](https://settylab.github.io/convert2anndata/reference/attach_reductions_seurat.md)
  : Attach obsm embeddings to a Seurat object as dimensional reductions
- [`attach_alt_experiments_sce()`](https://settylab.github.io/convert2anndata/reference/attach_alt_experiments_sce.md)
  : Attach altExps to SingleCellExperiment Object
- [`convert_graph_to_colPair()`](https://settylab.github.io/convert2anndata/reference/convert_graph_to_colPair.md)
  : Convert Seurat Graph to colPair-Compatible Format

## Internal building blocks

Exported for advanced use and testing, but not part of the recommended
user-facing API. Most users only need the converters in the first
section.

- [`align_metadata_seurat()`](https://settylab.github.io/convert2anndata/reference/align_metadata_seurat.md)
  : Align Metadata to Cell Names in Seurat Object
- [`convert_commands()`](https://settylab.github.io/convert2anndata/reference/convert_commands.md)
  : Convert SeuratCommand objects into a serializable format
- [`create_combined_assay()`](https://settylab.github.io/convert2anndata/reference/create_combined_assay.md)
  : Create Combined Assay in Seurat Object
- [`default_reduction_map()`](https://settylab.github.io/convert2anndata/reference/default_reduction_map.md)
  : Default mapping from AnnData obsm keys to Seurat reduction (name,
  key) pairs
- [`ensure_csparse_matrix()`](https://settylab.github.io/convert2anndata/reference/ensure_csparse_matrix.md)
  : Ensure Sparse Matrix is in a Format Compatible with AnnData
- [`extract_counts_matrix()`](https://settylab.github.io/convert2anndata/reference/extract_counts_matrix.md)
  : Extract Counts Matrix from Seurat Assay
- [`extract_data()`](https://settylab.github.io/convert2anndata/reference/extract_data.md)
  : Extract Data from SingleCellExperiment
- [`extract_pairs()`](https://settylab.github.io/convert2anndata/reference/extract_pairs.md)
  : Extract Pairwise Data from SingleCellExperiment
- [`identify_alt_exps_seurat()`](https://settylab.github.io/convert2anndata/reference/identify_alt_exps_seurat.md)
  : Identify Assays to Move to altExps in Seurat Object
- [`preserve_reductions_seurat()`](https://settylab.github.io/convert2anndata/reference/preserve_reductions_seurat.md)
  : Preserve Dimensional Reductions from Seurat Object
- [`process_alt_experiments()`](https://settylab.github.io/convert2anndata/reference/process_alt_experiments.md)
  : Process Alternative Experiments
- [`process_dimensional_reductions()`](https://settylab.github.io/convert2anndata/reference/process_dimensional_reductions.md)
  : Process Dimensional Reductions
- [`process_layers()`](https://settylab.github.io/convert2anndata/reference/process_layers.md)
  : Process Layers in Seurat Assay
- [`process_main_assay()`](https://settylab.github.io/convert2anndata/reference/process_main_assay.md)
  : Process Main Assay
- [`process_metadata_and_pairwise()`](https://settylab.github.io/convert2anndata/reference/process_metadata_and_pairwise.md)
  : Process Metadata and Pairwise Matrices
- [`process_other_assays()`](https://settylab.github.io/convert2anndata/reference/process_other_assays.md)
  : Process Other Assays
- [`timestamped_cat()`](https://settylab.github.io/convert2anndata/reference/timestamped_cat.md)
  : Print Messages with Timestamp
- [`update_seurat_object()`](https://settylab.github.io/convert2anndata/reference/update_seurat_object.md)
  : Update Seurat Object if Necessary
