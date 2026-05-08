#' Convert an AnnData object (or a .h5ad file) to a Seurat object
#'
#' Reverse direction of `convert_to_anndata()`: takes an AnnData object
#' (typically loaded with `anndata::read_h5ad()`) or a path to an `.h5ad`
#' file and constructs a Seurat object. Counts go into the `counts` layer,
#' `adata$X` into the `data` layer, `adata$obs` into the meta.data, and
#' entries of `adata$obsm` are attached as `DimReduc` objects.
#'
#' @param adata An AnnData object, OR a character path to an `.h5ad` file
#'   (in which case it is read with `anndata::read_h5ad()`).
#' @param counts_layer Name(s) of the layer in `adata$layers` to use as the
#'   counts matrix. May be a single string or a character vector of
#'   candidate names; the first one present in `adata$layers` is used.
#'   If none are present, `adata$X` is used as counts and no separate
#'   `data` layer is added. Defaults to
#'   `c("counts", "raw_counts", "raw_count")`.
#' @param assay Name of the resulting Seurat assay. Defaults to `"RNA"`.
#' @param reduction_map Optional named list passed to
#'   `attach_reductions_seurat()` to override or extend the default obsm-key
#'   to Seurat-reduction mapping. See `default_reduction_map()`.
#' @param conda_env Optional conda environment to activate before reading
#'   the file. Equivalent to calling `setup_anndata_python(conda_env)`
#'   first. Only relevant when `adata` is a path.
#' @param orig.ident Optional value for the `orig.ident` column of the
#'   resulting Seurat meta.data. If `NULL` (default), the value is taken
#'   from `adata$uns$conversion_source` if present, otherwise from the file
#'   basename when `adata` was passed as a path, otherwise `"AnnData"`.
#' @param use_raw How to handle `adata$raw` (a separate AnnData snapshot
#'   often used by scanpy to retain raw counts after normalising `X`).
#'   `"auto"` (default) uses `raw$X` as the counts layer when no
#'   `counts_layer` candidate matched and `adata$raw` is set; `TRUE`
#'   always uses `raw$X` as counts (and falls back to candidates / `X`
#'   if it is missing); `FALSE` ignores raw entirely.
#' @param attach_obsp Logical. When `TRUE` (default), entries of
#'   `adata$obsp` (cell x cell pairwise matrices, e.g. nearest-neighbor
#'   graphs from `sc.pp.neighbors`) are attached as Seurat `Graphs`.
#' @return A Seurat object.
#' @details The returned object has the AnnData orientation transposed: Seurat
#'   stores genes x cells while AnnData stores cells x genes. Cell and feature
#'   names are taken from `adata$obs_names` and `adata$var_names`.
#'
#'   `adata$obs` is attached as `meta.data`. `adata$var` is attached via
#'   `seurat_obj[["RNA"]][[]] <- ...` when present and dimensionally
#'   compatible.
#' @examples
#' \dontrun{
#' # Path entrypoint:
#' seurat_obj <- convert_anndata_to_seurat("pbmc.h5ad")
#'
#' # Object entrypoint:
#' adata <- anndata::read_h5ad("pbmc.h5ad")
#' seurat_obj <- convert_anndata_to_seurat(adata)
#'
#' # Custom layer + extra reductions:
#' seurat_obj <- convert_anndata_to_seurat(
#'   "pbmc.h5ad",
#'   counts_layer = c("counts", "raw_counts"),
#'   reduction_map = list(X_my_emb = list(name = "myemb", key = "MyEmb_"))
#' )
#' }
#' @importFrom Seurat CreateSeuratObject SetAssayData
#' @importFrom Matrix t
#' @importFrom methods as
#' @export
convert_anndata_to_seurat <- function(adata,
                                      counts_layer = c("counts", "raw_counts", "raw_count"),
                                      assay = "RNA",
                                      reduction_map = list(),
                                      conda_env = NULL,
                                      orig.ident = NULL,
                                      use_raw = c("auto", "always", "never"),
                                      attach_obsp = TRUE) {
  if (is.logical(use_raw)) {
    use_raw <- if (isTRUE(use_raw)) "always" else "never"
  }
  use_raw <- match.arg(use_raw)
  if (is.null(adata)) stop("`adata` is NULL.")

  source_path <- NULL
  if (is.character(adata) && length(adata) == 1) {
    source_path <- adata
    if (!file.exists(source_path)) {
      stop("Input file does not exist: ", source_path)
    }
    if (!is.null(conda_env)) setup_anndata_python(conda_env)
    # Fail fast with an actionable message if Python / anndata / numpy
    # aren't wired up correctly -- otherwise read_h5ad surfaces a
    # cryptic error from deep inside Python.
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

  raw_mat <- NULL
  used_raw <- FALSE
  if (use_raw %in% c("auto", "always")) {
    use_when_no_layer <- !has_counts_layer && use_raw == "auto"
    if (use_raw == "always" || use_when_no_layer) {
      raw_mat <- extract_anndata_raw(adata, obs_names = obs_names)
    }
  }

  if (!is.null(raw_mat)) {
    used_raw <- TRUE
    counts_mat <- raw_mat
    timestamped_cat(sprintf(
      "Using adata$raw$X (%d x %d) as counts.\n",
      nrow(raw_mat), ncol(raw_mat)
    ))
  } else {
    counts_mat <- extract_anndata_X(
      adata,
      layer = selected_layer,
      obs_names = obs_names,
      var_names = var_names
    )
  }

  meta.data <- NULL
  if (!is.null(adata$obs)) {
    meta.data <- tryCatch(as.data.frame(adata$obs), error = function(e) NULL)
    if (!is.null(meta.data)) {
      if (!"orig.ident" %in% colnames(meta.data)) {
        meta.data$orig.ident <- resolve_orig_ident(orig.ident, adata, source_path)
      }
      if (!is.null(obs_names)) rownames(meta.data) <- obs_names
    }
  }

  seurat_obj <- Seurat::CreateSeuratObject(
    counts = Matrix::t(counts_mat),
    meta.data = meta.data,
    assay = assay
  )

  # Add adata$X as the data layer when:
  #   - we used a counts layer for counts (X is the normalised companion), OR
  #   - we used raw and the main X has the same gene set as raw
  # If we used raw with a *different* gene set (typical scanpy pattern: raw
  # is unfiltered), attaching X would fail because the dims don't match the
  # Seurat object's gene axis. We skip in that case and emit a notice.
  if ((has_counts_layer || used_raw) && !is.null(adata$X)) {
    data_mat <- extract_anndata_X(adata, layer = NULL, obs_names = obs_names, var_names = var_names)
    data_mat <- Matrix::t(data_mat)
    if (nrow(data_mat) == nrow(seurat_obj) && ncol(data_mat) == ncol(seurat_obj)) {
      # Seurat sanitizes feature names (e.g. 'gene_01' -> 'gene-01') when the
      # object is created, so the data matrix's rownames/colnames must be
      # realigned to whatever the seurat object actually carries before
      # SetAssayData -- positions are preserved, only labels were rewritten.
      rownames(data_mat) <- rownames(seurat_obj)
      colnames(data_mat) <- colnames(seurat_obj)
      seurat_obj <- Seurat::SetAssayData(
        object = seurat_obj,
        assay = assay,
        layer = "data",
        new.data = data_mat
      )
      timestamped_cat(sprintf("Added adata$X to assay '%s' as 'data' layer.\n", assay))
    } else if (used_raw) {
      timestamped_cat(sprintf(
        "adata$X (%d x %d) has a different shape than adata$raw$X (%d x %d); skipping data layer.\n",
        nrow(data_mat), ncol(data_mat), nrow(seurat_obj), ncol(seurat_obj)
      ))
    }
  }

  obsm <- extract_anndata_obsm(adata, obs_names = obs_names)
  seurat_obj <- attach_reductions_seurat(
    seurat_obj, obsm, assay = assay, reduction_map = reduction_map
  )

  if (isTRUE(attach_obsp)) {
    seurat_obj <- attach_obsp_graphs(seurat_obj, adata, obs_names = obs_names,
                                     assay = assay)
  }

  # Choose the right `var` source. When we used raw, the active feature
  # space is raw.var; otherwise it's adata.var.
  var_source <- if (used_raw) {
    tryCatch(adata$raw$var, error = function(e) NULL)
  } else {
    tryCatch(adata$var, error = function(e) NULL)
  }
  if (!is.null(var_source)) {
    var_df <- tryCatch(as.data.frame(var_source), error = function(e) NULL)
    if (!is.null(var_df) && nrow(var_df) == nrow(seurat_obj)) {
      # Same realignment story as the data matrix: Seurat may have
      # sanitized feature names, and feature meta has to match exactly.
      rownames(var_df) <- rownames(seurat_obj)
      tryCatch(
        seurat_obj[[assay]][[]] <- var_df,
        error = function(e) {
          timestamped_cat("WARNING: could not attach var as feature meta:", conditionMessage(e), "\n")
        }
      )
    }
  }

  timestamped_cat("Conversion to Seurat complete.\n")
  seurat_obj
}

#' @keywords internal
attach_obsp_graphs <- function(seurat_obj, adata, obs_names = NULL, assay = "RNA") {
  obsp <- tryCatch(extract_anndata_obsp(adata, obs_names = obs_names),
                   error = function(e) list())
  if (length(obsp) == 0) return(seurat_obj)

  for (k in names(obsp)) {
    g <- obsp[[k]]
    if (nrow(g) != ncol(seurat_obj)) {
      warning(sprintf(
        "obsp['%s'] has %d cells but the Seurat object has %d; skipping.",
        k, nrow(g), ncol(seurat_obj)
      ))
      next
    }
    rownames(g) <- colnames(seurat_obj)
    colnames(g) <- colnames(seurat_obj)
    graph_obj <- tryCatch(
      SeuratObject::as.Graph(g),
      error = function(e) NULL
    )
    if (is.null(graph_obj)) {
      warning(sprintf("Could not coerce obsp['%s'] to a Seurat Graph; skipping.", k))
      next
    }
    SeuratObject::DefaultAssay(graph_obj) <- assay
    # Seurat convention is "<assay>_<name>" for graphs, so stay consistent.
    graph_name <- paste0(assay, "_", k)
    seurat_obj[[graph_name]] <- graph_obj
    timestamped_cat(sprintf("Attached obsp['%s'] as Seurat graph '%s'.\n", k, graph_name))
  }
  seurat_obj
}

pick_first_layer <- function(adata, candidates) {
  if (is.null(candidates)) return(NULL)
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) == 0) return(NULL)
  layer_keys <- anndata_mapping_keys(adata$layers)
  hit <- candidates[candidates %in% layer_keys]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

resolve_orig_ident <- function(orig.ident, adata, source_path) {
  if (!is.null(orig.ident)) return(orig.ident)
  src <- tryCatch(adata$uns[["conversion_source"]], error = function(e) NULL)
  if (!is.null(src) && length(src) == 1) return(as.character(src))
  if (!is.null(source_path)) {
    return(basename(tools::file_path_sans_ext(source_path)))
  }
  "AnnData"
}
