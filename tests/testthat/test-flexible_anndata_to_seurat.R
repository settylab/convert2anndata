library(testthat)
library(convert2anndata)
library(Matrix)
library(anndata)

skip_if_no_anndata <- function() {
  skip_if_not_installed("anndata")
  skip_if_not_installed("reticulate")
  if (!reticulate::py_available(initialize = TRUE)) {
    skip("Python not available for reticulate.")
  }
  if (!reticulate::py_module_available("anndata")) {
    skip("Python 'anndata' module not installed.")
  }
}

build_adata <- function(layer_names = "counts",
                       obsm_keys = c("X_pca"),
                       n_cells = 12, n_genes = 8) {
  X <- matrix(rpois(n_cells * n_genes, lambda = 3),
              nrow = n_cells, ncol = n_genes)
  rownames(X) <- paste0("cell", seq_len(n_cells))
  colnames(X) <- paste0("gene", seq_len(n_genes))

  layers <- if (length(layer_names)) {
    setNames(replicate(length(layer_names), X, simplify = FALSE), layer_names)
  } else NULL

  obsm <- if (length(obsm_keys)) {
    setNames(
      lapply(obsm_keys, function(k) matrix(rnorm(n_cells * 3), n_cells, 3)),
      obsm_keys
    )
  } else NULL

  AnnData(X = X, obs = data.frame(cell_type = rep(c("A", "B"), n_cells / 2),
                                  row.names = rownames(X)),
          var = data.frame(gene_type = rep("g", n_genes), row.names = colnames(X)),
          layers = layers, obsm = obsm)
}

test_that("counts_layer accepts a vector of candidates", {
  skip_if_no_anndata()
  skip_if_not_installed("Seurat")
  ad <- build_adata(layer_names = "raw_counts")
  s <- convert_anndata_to_seurat(
    ad,
    counts_layer = c("counts", "raw_counts", "raw_count")
  )
  expect_s4_class(s, "Seurat")
  counts <- as.matrix(Seurat::GetAssayData(s, layer = "counts"))
  raw <- as.matrix(ad$layers[["raw_counts"]])
  expect_equal(unname(counts), unname(t(raw)))
})

test_that("counts_layer falls back to X when no candidate matches", {
  skip_if_no_anndata()
  skip_if_not_installed("Seurat")
  ad <- build_adata(layer_names = "logcounts")  # neither counts nor raw_counts
  s <- convert_anndata_to_seurat(
    ad, counts_layer = c("counts", "raw_counts")
  )
  expect_s4_class(s, "Seurat")
  # No data layer should be added because we used X as counts.
  layers <- SeuratObject::Layers(s[["RNA"]])
  expect_false("data" %in% layers && length(layers) > 1 &&
               !identical(unname(as.matrix(Seurat::GetAssayData(s, layer = "data"))),
                          unname(as.matrix(Seurat::GetAssayData(s, layer = "counts")))))
})

test_that("convert_anndata_to_seurat accepts a file path and reads it", {
  skip_if_no_anndata()
  skip_if_not_installed("Seurat")
  ad <- build_adata()
  path <- tempfile(fileext = ".h5ad")
  write_h5ad(ad, path)

  s <- convert_anndata_to_seurat(path, counts_layer = "counts")
  expect_s4_class(s, "Seurat")
  expect_equal(ncol(s), 12)
  expect_equal(nrow(s), 8)
  # orig.ident should be the file basename when no uns$conversion_source is set
  expect_true("orig.ident" %in% colnames(s@meta.data))
  expect_equal(unique(as.character(s@meta.data$orig.ident)),
               tools::file_path_sans_ext(basename(path)))
})

test_that("convert_anndata_to_seurat respects an explicit orig.ident", {
  skip_if_no_anndata()
  skip_if_not_installed("Seurat")
  ad <- build_adata()
  s <- convert_anndata_to_seurat(ad, counts_layer = "counts",
                                 orig.ident = "my_sample")
  expect_equal(unique(as.character(s@meta.data$orig.ident)), "my_sample")
})

test_that("missing input file raises a clear error", {
  expect_error(
    convert_anndata_to_seurat("/no/such/file.h5ad"),
    "does not exist"
  )
})

test_that("default_reduction_map covers common scanpy embeddings", {
  m <- default_reduction_map()
  expect_true(all(c("X_pca", "X_umap", "X_tsne") %in% names(m)))
  expect_equal(m$X_pca$name, "pca")
  expect_equal(m$X_pca$key, "PC_")
})

test_that("reduction_map overrides defaults and adds new keys", {
  skip_if_no_anndata()
  skip_if_not_installed("Seurat")
  ad <- build_adata(obsm_keys = c("X_pca", "X_my_emb"))
  custom <- list(
    X_pca    = list(name = "PCA_renamed", key = "Renamed_"),
    X_my_emb = list(name = "myemb",       key = "MyEmb_")
  )
  s <- convert_anndata_to_seurat(ad, counts_layer = "counts",
                                 reduction_map = custom)
  expect_true("PCA_renamed" %in% names(s@reductions))
  expect_true("myemb" %in% names(s@reductions))
  expect_equal(s@reductions$PCA_renamed@key, "Renamed_")
  expect_equal(s@reductions$myemb@key, "MyEmb_")
})

test_that("reduction_map can disable a default mapping with name = NA", {
  skip_if_no_anndata()
  skip_if_not_installed("Seurat")
  ad <- build_adata(obsm_keys = c("X_pca", "X_umap"))
  s <- convert_anndata_to_seurat(
    ad, counts_layer = "counts",
    reduction_map = list(X_umap = list(name = NA, key = NA))
  )
  expect_true("pca" %in% names(s@reductions))
  expect_false("umap" %in% names(s@reductions))
})

test_that("convert_anndata_to_sce accepts a file path", {
  skip_if_no_anndata()
  ad <- build_adata()
  path <- tempfile(fileext = ".h5ad")
  write_h5ad(ad, path)

  sce <- convert_anndata_to_sce(path, counts_layer = "counts")
  expect_s4_class(sce, "SingleCellExperiment")
})

test_that("convert_anndata_to_sce counts_layer accepts a vector of candidates", {
  skip_if_no_anndata()
  ad <- build_adata(layer_names = "raw_count")
  sce <- convert_anndata_to_sce(
    ad, counts_layer = c("counts", "raw_count")
  )
  expect_s4_class(sce, "SingleCellExperiment")
  expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
})

test_that("pick_first_layer returns NULL on empty/NULL candidates", {
  skip_if_no_anndata()
  ad <- build_adata()
  expect_null(convert2anndata:::pick_first_layer(ad, NULL))
  expect_null(convert2anndata:::pick_first_layer(ad, character(0)))
  expect_null(convert2anndata:::pick_first_layer(ad, ""))
})

test_that("orig.ident derives from uns$conversion_source when present", {
  skip_if_no_anndata()
  skip_if_not_installed("Seurat")
  ad <- build_adata()
  ad$uns[["conversion_source"]] <- "labeled_in_uns"
  s <- convert_anndata_to_seurat(ad, counts_layer = "counts")
  expect_equal(unique(as.character(s@meta.data$orig.ident)), "labeled_in_uns")
})
