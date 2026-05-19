#' Default mapping from AnnData obsm keys to Seurat reduction (name, key) pairs
#'
#' Provides the canonical Seurat naming for the obsm keys conventionally
#' produced by scanpy. Users who want extra mappings should pass a list with
#' the same shape into `attach_reductions_seurat()` (or its callers).
#'
#' @return A named list. Each name is an obsm key (e.g. `X_pca`); each
#'   element is a list with `name` (Seurat reduction name, e.g. `"pca"`)
#'   and `key` (column-name prefix, e.g. `"PC_"`).
#' @export
default_reduction_map <- function() {
  list(
    X_pca       = list(name = "pca",       key = "PC_"),
    X_umap      = list(name = "umap",      key = "UMAP_"),
    X_tsne      = list(name = "tsne",      key = "tSNE_"),
    X_diffmap   = list(name = "diffmap",   key = "DM_"),
    X_phate     = list(name = "phate",     key = "PHATE_"),
    X_harmony   = list(name = "harmony",   key = "harmony_"),
    X_scvi      = list(name = "scvi",      key = "scVI_"),
    X_lsi       = list(name = "lsi",       key = "LSI_"),
    X_spectral  = list(name = "spectral",  key = "spectral_")
  )
}

#' Attach obsm embeddings to a Seurat object as dimensional reductions
#'
#' Uses an obsm-key -> (reduction-name, key-prefix) mapping. Keys not in the
#' mapping fall back to a derived name: leading `X_` stripped, lowercased,
#' with a sanitized key prefix.
#'
#' @param seurat_obj A Seurat object to attach reductions to.
#' @param obsm Named list of numeric matrices (cells x dim), as returned by
#'   `extract_anndata_obsm()`.
#' @param assay Assay name to associate the reductions with. Defaults to "RNA".
#' @param reduction_map Named list overriding or extending
#'   `default_reduction_map()`. Each entry maps an obsm key (e.g. `"X_pca"`)
#'   to a list with `name` and `key`. The supplied map is merged on top of
#'   the defaults; pass an empty list to keep defaults, or pass an entry
#'   with `name = NA` to disable a default.
#' @return The Seurat object with reductions attached.
#' @importFrom Seurat CreateDimReducObject
#' @export
attach_reductions_seurat <- function(seurat_obj, obsm, assay = "RNA",
                                     reduction_map = list()) {
  if (length(obsm) == 0) return(seurat_obj)

  merged_map <- default_reduction_map()
  for (k in names(reduction_map)) {
    merged_map[[k]] <- reduction_map[[k]]
  }

  for (k in names(obsm)) {
    mat <- obsm[[k]]
    if (nrow(mat) != ncol(seurat_obj)) {
      warning(sprintf(
        "obsm['%s'] has %d rows but the Seurat object has %d cells; skipping.",
        k, nrow(mat), ncol(seurat_obj)
      ))
      next
    }
    if (k %in% names(merged_map)) {
      entry <- merged_map[[k]]
      if (is.null(entry) || (length(entry$name) == 1 && is.na(entry$name))) {
        timestamped_cat(sprintf("Skipping obsm['%s'] (disabled in reduction_map).\n", k))
        next
      }
      red_name <- entry$name
      red_key  <- entry$key %||% paste0(red_name, "_")
    } else {
      red_name <- tolower(sub("^X_", "", k))
      red_key  <- paste0(red_name, "_")
    }
    sanitized_key <- gsub("[^[:alnum:]]", "", red_key)
    if (!grepl("_$", sanitized_key)) sanitized_key <- paste0(sanitized_key, "_")
    colnames(mat) <- paste0(sanitized_key, seq_len(ncol(mat)))
    # Seurat may have rewritten cell names during object construction; align
    # the embedding rownames to the seurat object so CreateDimReducObject
    # accepts them.
    rownames(mat) <- colnames(seurat_obj)

    seurat_obj[[red_name]] <- Seurat::CreateDimReducObject(
      embeddings = mat,
      key = sanitized_key,
      assay = assay
    )
    timestamped_cat(sprintf("Added obsm['%s'] as reduction '%s'.\n", k, red_name))
  }
  seurat_obj
}

`%||%` <- function(a, b) if (is.null(a)) b else a
