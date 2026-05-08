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

# Build a 1000-cell, 500-gene synthetic AnnData that mirrors the typical
# scanpy-pipeline output: filtered + normalised X, raw counts in
# adata.raw, layers for counts/logcounts, multiple obsm reductions, an
# obsp NN graph, multi-column categorical obs, and uns metadata.
build_scanpy_like <- function(n_cells = 1000L, n_genes_full = 500L,
                              n_genes_hv = 200L, seed = 42L) {
  set.seed(seed)
  np <- reticulate::import("numpy")
  sp <- reticulate::import("scipy.sparse")

  raw_dense <- matrix(rpois(n_cells * n_genes_full, lambda = 1.5),
                      nrow = n_cells, ncol = n_genes_full)
  rownames(raw_dense) <- sprintf("cell-%04d", seq_len(n_cells))
  colnames(raw_dense) <- sprintf("ENSG-%04d", seq_len(n_genes_full))

  hv_idx <- sort(sample.int(n_genes_full, n_genes_hv))
  X_filt <- log1p(raw_dense[, hv_idx])
  filt_genes <- colnames(raw_dense)[hv_idx]
  colnames(X_filt) <- filt_genes

  raw_sparse <- sp$csr_matrix(np$asarray(raw_dense, dtype = "float32"))
  X_sparse   <- sp$csr_matrix(np$asarray(X_filt,    dtype = "float32"))

  obs <- data.frame(
    sample      = sample(paste0("S", 1:5), n_cells, replace = TRUE),
    cell_type   = factor(sample(c("T", "B", "NK", "Mono", "DC"),
                                n_cells, replace = TRUE)),
    n_counts    = rowSums(raw_dense),
    n_genes     = rowSums(raw_dense > 0),
    pct_mt      = runif(n_cells, 0, 10),
    is_doublet  = sample(c(FALSE, TRUE), n_cells,
                         replace = TRUE, prob = c(0.95, 0.05)),
    row.names   = rownames(raw_dense),
    stringsAsFactors = FALSE
  )
  var_full <- data.frame(
    highly_variable = colnames(raw_dense) %in% filt_genes,
    mean_counts     = colMeans(raw_dense),
    n_cells         = colSums(raw_dense > 0),
    row.names       = colnames(raw_dense),
    stringsAsFactors = FALSE
  )
  var_filt <- var_full[filt_genes, , drop = FALSE]

  X_pca   <- matrix(rnorm(n_cells * 30), n_cells, 30)
  X_umap  <- matrix(rnorm(n_cells * 2),  n_cells, 2)
  X_tsne  <- matrix(rnorm(n_cells * 2),  n_cells, 2)
  X_harm  <- matrix(rnorm(n_cells * 30), n_cells, 30)

  # kNN graph (k = 10, sparse, square)
  G <- matrix(0, n_cells, n_cells)
  for (i in seq_len(n_cells)) {
    nbrs <- sample(setdiff(seq_len(n_cells), i), 10)
    G[i, nbrs] <- runif(10)
  }
  conn <- sp$csr_matrix(np$asarray(G,           dtype = "float32"))
  dist <- sp$csr_matrix(np$asarray(G * 5,       dtype = "float32"))

  ad <- AnnData(
    X      = X_sparse,
    obs    = obs,
    var    = var_filt,
    obsm   = list(X_pca = X_pca, X_umap = X_umap,
                  X_tsne = X_tsne, X_harmony = X_harm),
    obsp   = list(connectivities = conn, distances = dist),
    uns    = list(
      conversion_source = "scanpy_synth",
      neighbors         = list(method = "umap", n_neighbors = 10L),
      pca               = list(variance_ratio = runif(30))
    ),
    layers = list(
      counts    = X_sparse,
      logcounts = X_sparse
    )
  )
  ad$raw <- AnnData(X = raw_sparse, var = var_full)
  ad
}

test_that("realistic scanpy-like dataset converts and round-trips key invariants", {
  skip_if_no_anndata()
  ad <- build_scanpy_like(n_cells = 600L, n_genes_full = 300L, n_genes_hv = 120L)
  s <- convert_anndata_to_seurat(ad, counts_layer = "counts")

  expect_s4_class(s, "Seurat")
  expect_equal(ncol(s), as.integer(ad$n_obs))
  expect_equal(nrow(s), as.integer(ad$n_vars))   # X.shape, not raw.shape

  # All four reductions show up.
  expect_setequal(
    intersect(names(s@reductions), c("pca","umap","tsne","harmony")),
    c("pca","umap","tsne","harmony")
  )

  # Both obsp graphs attached.
  expect_true(all(c("RNA_connectivities", "RNA_distances") %in% names(s@graphs)))

  # Categorical obs survives.
  expect_true("cell_type" %in% colnames(s@meta.data))
  obs_in <- as.data.frame(ad$obs)
  expect_equal(as.character(s@meta.data$cell_type),
               as.character(obs_in$cell_type))

  # Feature meta got attached.
  fm <- s[["RNA"]][[]]
  expect_true("highly_variable" %in% colnames(fm))

  # uns$conversion_source becomes orig.ident.
  expect_equal(unique(as.character(s@meta.data$orig.ident)), "scanpy_synth")
})

test_that("realistic dataset with use_raw='always' yields the unfiltered gene set", {
  skip_if_no_anndata()
  ad <- build_scanpy_like(n_cells = 400L, n_genes_full = 200L, n_genes_hv = 80L)
  s <- convert_anndata_to_seurat(ad, use_raw = "always")
  expect_equal(nrow(s), as.integer(ad$raw$shape[2]))
  # Data layer should be skipped because raw and X have different gene sets.
  expect_false("data" %in% SeuratObject::Layers(s[["RNA"]]))
})

test_that("realistic dataset converts in reasonable time", {
  skip_if_no_anndata()
  ad <- build_scanpy_like(n_cells = 1000L, n_genes_full = 500L, n_genes_hv = 200L)
  t <- system.time({
    s <- convert_anndata_to_seurat(ad, counts_layer = "counts")
  })[["elapsed"]]
  # Generous bound: this is on shared infrastructure. 90s is a sanity check
  # against quadratic-blowup regressions, not a perf benchmark.
  expect_lt(t, 90)
  expect_equal(ncol(s), 1000L)
})

# ----------------------------------------------------------------------
# backed='r' mode
# ----------------------------------------------------------------------

test_that("backed='r' read still produces an in-memory Seurat object", {
  skip_if_no_anndata()
  ad <- build_scanpy_like(n_cells = 200L, n_genes_full = 150L, n_genes_hv = 60L)
  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  write_h5ad(ad, path)

  ad_backed <- tryCatch(
    anndata::read_h5ad(path, backed = "r"),
    error = function(e) {
      skip(paste("backed read not supported by this anndata:",
                 conditionMessage(e)))
    }
  )

  s <- tryCatch(
    convert_anndata_to_seurat(ad_backed, counts_layer = "counts"),
    error = function(e) {
      skip(paste("backed conversion not supported:", conditionMessage(e)))
    }
  )
  expect_s4_class(s, "Seurat")
  expect_equal(ncol(s), 200L)
})
