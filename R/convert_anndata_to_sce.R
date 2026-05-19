#' Convert an AnnData object (or a .h5ad file) to a SingleCellExperiment
#'
#' Reverse direction of `convert_to_anndata()` for the SCE branch.
#' `adata$X` (or `adata$layers[[counts_layer]]` if present) becomes the
#' `counts` assay, additional layers become other assays, `adata$obs` is
#' written to `colData`, `adata$var` to `rowData`, and `adata$obsm` entries
#' to `reducedDims`.
#'
#' @param adata An AnnData object, OR a character path to an `.h5ad` file
#'   (in which case it is read with `anndata::read_h5ad()`).
#' @param counts_layer Name(s) of the layer in `adata$layers` to use as the
#'   `counts` assay. May be a single string or a character vector of
#'   candidates; the first one present in `adata$layers` is used. If none
#'   are present, `adata$X` is used as `counts`. Defaults to
#'   `c("counts", "raw_counts", "raw_count")`.
#' @param conda_env Optional conda environment to activate before reading
#'   the file. Equivalent to calling `setup_anndata_python(conda_env)`
#'   first. Only relevant when `adata` is a path.
#' @return A SingleCellExperiment object.
#' @examples
#' \dontrun{
#' # Path entrypoint:
#' sce <- convert_anndata_to_sce("pbmc.h5ad")
#' # Object entrypoint:
#' adata <- anndata::read_h5ad("pbmc.h5ad")
#' sce <- convert_anndata_to_sce(adata)
#' }
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims<-
#' @importFrom SummarizedExperiment colData<- rowData<-
#' @importFrom Matrix t
#' @importFrom methods as
#' @export
convert_anndata_to_sce <- function(adata,
                                   counts_layer = c("counts", "raw_counts", "raw_count"),
                                   conda_env = NULL) {
  if (is.null(adata)) stop("`adata` is NULL.")

  if (is.character(adata) && length(adata) == 1) {
    source_path <- adata
    if (!file.exists(source_path)) {
      stop("Input file does not exist: ", source_path)
    }
    if (!is.null(conda_env)) setup_anndata_python(conda_env)
    check_anndata_python()
    timestamped_cat(sprintf("Reading AnnData from '%s'.\n", source_path))
    adata <- tryCatch(
      anndata::read_h5ad(source_path),
      error = function(e) stop("Failed to read AnnData h5ad file: ", conditionMessage(e))
    )
  }

  timestamped_cat("Summary of AnnData object:\n\n")
  print(adata)
  cat("\n")

  names <- anndata_names(adata)
  obs_names <- names$obs_names
  var_names <- names$var_names

  selected_layer <- pick_first_layer(adata, counts_layer)
  has_counts_layer <- !is.null(selected_layer)

  counts_mat <- extract_anndata_X(
    adata,
    layer = selected_layer,
    obs_names = obs_names,
    var_names = var_names
  )
  # SCE/Seurat orientation: features x cells; AnnData is cells x features.
  counts_mat <- Matrix::t(counts_mat)

  assays_list <- list(counts = counts_mat)

  if (has_counts_layer && !is.null(adata$X)) {
    X_mat <- extract_anndata_X(adata, layer = NULL, obs_names = obs_names, var_names = var_names)
    assays_list[["X"]] <- Matrix::t(X_mat)
  }

  exclude_layers <- if (has_counts_layer) selected_layer else character(0)
  other_layers <- extract_anndata_layers(
    adata,
    obs_names = obs_names,
    var_names = var_names,
    exclude = exclude_layers
  )
  for (k in names(other_layers)) {
    assays_list[[k]] <- Matrix::t(other_layers[[k]])
  }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = assays_list)

  if (!is.null(adata$obs)) {
    obs_df <- tryCatch(as.data.frame(adata$obs), error = function(e) NULL)
    if (!is.null(obs_df) && nrow(obs_df) == ncol(sce)) {
      if (!is.null(obs_names)) rownames(obs_df) <- obs_names
      SummarizedExperiment::colData(sce) <- methods::as(obs_df, "DataFrame")
    }
  }
  if (!is.null(adata$var)) {
    var_df <- tryCatch(as.data.frame(adata$var), error = function(e) NULL)
    if (!is.null(var_df) && nrow(var_df) == nrow(sce)) {
      if (!is.null(var_names)) rownames(var_df) <- var_names
      SummarizedExperiment::rowData(sce) <- methods::as(var_df, "DataFrame")
    }
  }

  obsm <- extract_anndata_obsm(adata, obs_names = obs_names)
  if (length(obsm) > 0) {
    # Strip leading "X_" so reducedDimNames are clean (e.g. "PCA", "UMAP").
    clean_names <- sub("^X_", "", names(obsm))
    names(obsm) <- toupper(clean_names)
    SingleCellExperiment::reducedDims(sce) <- obsm
  }

  timestamped_cat("Conversion to SingleCellExperiment complete.\n")
  sce
}
