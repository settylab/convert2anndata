library(testthat)
library(convert2anndata)
library(Matrix)
library(anndata)
library(Seurat)

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

build_X <- function(n_cells, n_genes) {
  X <- matrix(rpois(n_cells * n_genes, 3),
              nrow = n_cells, ncol = n_genes)
  rownames(X) <- sprintf("cell-%03d", seq_len(n_cells))
  colnames(X) <- sprintf("gene-%03d", seq_len(n_genes))
  X
}

# ----------------------------------------------------------------------
# adata.raw
# ----------------------------------------------------------------------

test_that("use_raw='auto' uses raw when no counts_layer matches", {
  skip_if_no_anndata()
  X <- build_X(12, 8)
  log_X <- log1p(X)
  ad <- AnnData(X = log_X)
  ad$raw <- AnnData(X = X)

  s <- convert_anndata_to_seurat(ad, counts_layer = c("counts","raw_counts"))
  counts <- as.matrix(GetAssayData(s, layer = "counts"))
  expect_equal(unname(counts), unname(t(X)))
  data <- as.matrix(GetAssayData(s, layer = "data"))
  expect_equal(unname(data), unname(t(log_X)))
})

test_that("use_raw='auto' does NOT override an existing counts_layer match", {
  skip_if_no_anndata()
  X <- build_X(10, 6)
  log_X <- log1p(X)
  layer_X <- X * 7L
  ad <- AnnData(X = log_X, layers = list(counts = layer_X))
  ad$raw <- AnnData(X = X)

  s <- convert_anndata_to_seurat(ad, counts_layer = "counts")
  counts <- as.matrix(GetAssayData(s, layer = "counts"))
  # The layer matched, so it wins -- raw should be ignored.
  expect_equal(unname(counts), unname(t(layer_X)))
})

test_that("use_raw='always' uses raw even when a counts_layer matches", {
  skip_if_no_anndata()
  X <- build_X(10, 6)
  log_X <- log1p(X)
  layer_X <- X * 7L
  ad <- AnnData(X = log_X, layers = list(counts = layer_X))
  ad$raw <- AnnData(X = X)

  s <- convert_anndata_to_seurat(ad, counts_layer = "counts", use_raw = "always")
  counts <- as.matrix(GetAssayData(s, layer = "counts"))
  expect_equal(unname(counts), unname(t(X)))
})

test_that("use_raw='never' ignores raw entirely", {
  skip_if_no_anndata()
  X <- build_X(10, 6)
  ad <- AnnData(X = log1p(X))
  ad$raw <- AnnData(X = X)

  s <- convert_anndata_to_seurat(ad, counts_layer = c("counts"), use_raw = "never")
  counts <- as.matrix(GetAssayData(s, layer = "counts"))
  # Raw was ignored; no counts layer matched; X (log1p) became counts.
  expect_equal(unname(counts), unname(t(log1p(X))))
})

test_that("use_raw can be passed as logical for back-compat", {
  skip_if_no_anndata()
  X <- build_X(10, 6)
  ad <- AnnData(X = log1p(X))
  ad$raw <- AnnData(X = X)

  s_t <- convert_anndata_to_seurat(ad, use_raw = TRUE)
  s_f <- convert_anndata_to_seurat(ad, use_raw = FALSE)
  expect_equal(unname(as.matrix(GetAssayData(s_t, layer = "counts"))), unname(t(X)))
  expect_equal(unname(as.matrix(GetAssayData(s_f, layer = "counts"))), unname(t(log1p(X))))
})

test_that("raw with a different gene set: counts uses raw shape; data layer skipped", {
  skip_if_no_anndata()
  # Typical scanpy pattern: adata is filtered to highly-variable genes
  # (here 6 of 12); adata.raw retains all 12 genes with raw counts.
  X_full <- build_X(10, 12)
  hv <- 1:6
  X_filt <- log1p(X_full[, hv])
  colnames(X_filt) <- colnames(X_full)[hv]

  ad <- AnnData(X = X_filt)
  ad$raw <- AnnData(X = X_full)

  s <- convert_anndata_to_seurat(ad, use_raw = "always")
  expect_equal(nrow(s), ncol(X_full))   # full gene set from raw
  expect_equal(ncol(s), nrow(X_full))   # cells unchanged
  # Data layer was skipped because the shapes differ.
  expect_false("data" %in% SeuratObject::Layers(s[["RNA"]]))
})

test_that("raw without an obs/var attached still produces a valid Seurat object", {
  skip_if_no_anndata()
  X <- build_X(8, 5)
  ad <- AnnData(X = log1p(X))
  ad$raw <- AnnData(X = X)
  s <- convert_anndata_to_seurat(ad, use_raw = "always")
  expect_s4_class(s, "Seurat")
  expect_equal(ncol(s), 8)
})

# ----------------------------------------------------------------------
# adata.obsp -> Seurat Graphs
# ----------------------------------------------------------------------

build_obsp_graph <- function(n) {
  np <- reticulate::import("numpy")
  sp <- reticulate::import("scipy.sparse")
  G <- matrix(0, n, n)
  set.seed(0)
  for (i in seq_len(n)) {
    nbrs <- sample(setdiff(seq_len(n), i), 2)
    G[i, nbrs] <- runif(2)
  }
  sp$csr_matrix(np$asarray(G, dtype = "float64"))
}

test_that("obsp connectivities and distances become Seurat Graphs", {
  skip_if_no_anndata()
  n <- 12
  X <- build_X(n, 8)
  ad <- AnnData(X = X)
  ad$obsp <- list(
    connectivities = build_obsp_graph(n),
    distances      = build_obsp_graph(n)
  )

  s <- convert_anndata_to_seurat(ad)
  graph_names <- names(s@graphs)
  expect_true("RNA_connectivities" %in% graph_names)
  expect_true("RNA_distances" %in% graph_names)
  g <- s@graphs[["RNA_connectivities"]]
  expect_equal(nrow(g), n)
  expect_equal(ncol(g), n)
  expect_equal(rownames(g), colnames(s))
})

test_that("obsp uses the custom assay name in graph names", {
  skip_if_no_anndata()
  n <- 8
  X <- build_X(n, 5)
  ad <- AnnData(X = X)
  ad$obsp <- list(connectivities = build_obsp_graph(n))
  s <- convert_anndata_to_seurat(ad, assay = "ATAC")
  expect_true("ATAC_connectivities" %in% names(s@graphs))
})

test_that("attach_obsp = FALSE skips obsp attachment", {
  skip_if_no_anndata()
  n <- 8
  X <- build_X(n, 5)
  ad <- AnnData(X = X)
  ad$obsp <- list(connectivities = build_obsp_graph(n))
  s <- convert_anndata_to_seurat(ad, attach_obsp = FALSE)
  expect_length(s@graphs, 0)
})

test_that("extract_anndata_obsp guards against non-square or shape-mismatched matrices", {
  skip_if_no_anndata()
  n <- 8
  # Use a fake adata-like environment that exposes obsp via its dollar-list.
  # We build a non-square matrix directly and route it through the helper.
  fake_obsp <- list(
    bad_rect   = Matrix::Matrix(matrix(runif(n * (n - 1)), n, n - 1), sparse = TRUE),
    good_square = Matrix::Matrix(matrix(runif(n * n), n, n), sparse = TRUE)
  )
  fake_adata <- list(obsp = fake_obsp)

  # Stub anndata_mapping_keys to read names() instead of py-keys()
  res <- expect_warning(
    extract_anndata_obsp(fake_adata, obs_names = sprintf("c%02d", seq_len(n))),
    regexp = "not square"
  )
  expect_named(res, "good_square")
  expect_equal(nrow(res$good_square), n)
})
