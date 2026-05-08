library(testthat)
library(convert2anndata)
library(Matrix)
library(anndata)
library(Seurat)

# Comprehensive coverage of AnnData -> Seurat conversion.
#
# Each test builds a small AnnData with one feature exercised in
# isolation, then asserts that the Seurat object preserves the
# semantics that a downstream analyst would expect.

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

# Repeatable RNG without polluting the user's global stream.
with_seed <- function(seed, expr) {
  old <- .Random.seed
  on.exit(assign(".Random.seed", old, envir = .GlobalEnv))
  set.seed(seed)
  force(expr)
}

build_dense_X <- function(n_cells, n_genes, lambda = 2.5, seed = 0L) {
  with_seed(seed, {
    X <- matrix(rpois(n_cells * n_genes, lambda),
                nrow = n_cells, ncol = n_genes)
  })
  # Names use hyphens, not underscores: Seurat rewrites underscores in
  # feature names to dashes (a documented Seurat policy), which would
  # confound assertions on label preservation. Cells are not sanitized
  # but we keep them consistent.
  rownames(X) <- sprintf("cell-%03d", seq_len(n_cells))
  colnames(X) <- sprintf("gene-%03d", seq_len(n_genes))
  X
}

# ----------------------------------------------------------------------
# A. Matrix variants: dense, CSR, CSC, dtypes
# ----------------------------------------------------------------------

test_that("dense X is preserved through the conversion", {
  skip_if_no_anndata()
  X <- build_dense_X(15, 10)
  ad <- AnnData(X = X)
  s <- convert_anndata_to_seurat(ad, counts_layer = "counts")
  counts <- as.matrix(GetAssayData(s, layer = "counts"))
  expect_equal(unname(counts), unname(t(X)))
  expect_equal(rownames(s), colnames(X))
  expect_equal(colnames(s), rownames(X))
})

test_that("CSR sparse X is preserved", {
  skip_if_no_anndata()
  sp <- reticulate::import("scipy.sparse")
  np <- reticulate::import("numpy")
  X <- build_dense_X(20, 12)
  X_py <- np$asarray(X, dtype = "float32")
  ad <- AnnData(X = sp$csr_matrix(X_py))
  s <- convert_anndata_to_seurat(ad)
  counts <- as.matrix(GetAssayData(s, layer = "counts"))
  expect_equal(unname(counts), unname(t(X)))
})

test_that("CSC sparse X is preserved", {
  skip_if_no_anndata()
  sp <- reticulate::import("scipy.sparse")
  np <- reticulate::import("numpy")
  X <- build_dense_X(20, 12)
  X_py <- np$asarray(X, dtype = "float32")
  ad <- AnnData(X = sp$csc_matrix(X_py))
  s <- convert_anndata_to_seurat(ad)
  counts <- as.matrix(GetAssayData(s, layer = "counts"))
  expect_equal(unname(counts), unname(t(X)))
})

test_that("float vs int dtype both round-trip", {
  skip_if_no_anndata()
  np <- reticulate::import("numpy")
  X <- build_dense_X(10, 6)
  ad_float <- AnnData(X = np$asarray(X, dtype = "float64"))
  ad_int   <- AnnData(X = np$asarray(X, dtype = "int32"))
  s_f <- convert_anndata_to_seurat(ad_float)
  s_i <- convert_anndata_to_seurat(ad_int)
  expect_equal(
    unname(as.matrix(GetAssayData(s_f, layer = "counts"))),
    unname(as.matrix(GetAssayData(s_i, layer = "counts")))
  )
})

# ----------------------------------------------------------------------
# B. Metadata fidelity (obs / var)
# ----------------------------------------------------------------------

test_that("obs preserves character, numeric, logical, and categorical columns", {
  skip_if_no_anndata()
  X <- build_dense_X(12, 8)
  obs <- data.frame(
    sample_id = paste0("S", seq_len(12)),
    n_counts  = rowSums(X),
    is_ctrl   = rep(c(TRUE, FALSE), 6),
    cell_type = factor(rep(c("A", "B", "C"), 4),
                       levels = c("A", "B", "C")),
    stringsAsFactors = FALSE,
    row.names = rownames(X)
  )
  ad <- AnnData(X = X, obs = obs)
  s <- convert_anndata_to_seurat(ad)
  md <- s@meta.data

  expect_setequal(
    setdiff(colnames(md), c("orig.ident", "nCount_RNA", "nFeature_RNA")),
    c("sample_id", "n_counts", "is_ctrl", "cell_type")
  )
  expect_equal(as.character(md$sample_id), as.character(obs$sample_id))
  expect_equal(unname(md$n_counts), unname(obs$n_counts))
  expect_equal(as.logical(md$is_ctrl), obs$is_ctrl)
  # cell_type round-trips through pandas Categorical; we only require
  # that the levels and per-cell labels are preserved, not the storage class.
  expect_equal(as.character(md$cell_type), as.character(obs$cell_type))
})

test_that("obs with NA values in numeric and categorical survive", {
  skip_if_no_anndata()
  X <- build_dense_X(12, 8)
  obs <- data.frame(
    n_counts = c(NA_real_, rowSums(X)[-1]),
    cell_type = factor(c(rep("A", 6), rep(NA, 6)), levels = c("A", "B")),
    row.names = rownames(X),
    stringsAsFactors = FALSE
  )
  ad <- AnnData(X = X, obs = obs)
  s <- convert_anndata_to_seurat(ad)
  md <- s@meta.data
  expect_true(is.na(md$n_counts[1]))
  expect_true(any(is.na(md$cell_type)))
})

test_that("an explicit orig.ident column in obs is preserved verbatim", {
  skip_if_no_anndata()
  X <- build_dense_X(8, 5)
  obs <- data.frame(
    orig.ident = paste0("donor_", rep(1:2, each = 4)),
    row.names  = rownames(X),
    stringsAsFactors = FALSE
  )
  ad <- AnnData(X = X, obs = obs)
  s <- convert_anndata_to_seurat(ad)
  expect_equal(as.character(s@meta.data$orig.ident), obs$orig.ident)
})

test_that("var attaches to the assay's feature meta and aligns to var_names", {
  skip_if_no_anndata()
  X <- build_dense_X(10, 6)
  var <- data.frame(
    gene_class = c("protein_coding", "lncRNA", "lncRNA", "protein_coding",
                   "protein_coding", "miRNA"),
    highly_variable = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
    row.names = colnames(X),
    stringsAsFactors = FALSE
  )
  ad <- AnnData(X = X, var = var)
  s <- convert_anndata_to_seurat(ad)
  fm <- s[["RNA"]][[]]
  expect_true("gene_class" %in% colnames(fm))
  expect_true("highly_variable" %in% colnames(fm))
  expect_equal(rownames(fm), colnames(X))
  expect_equal(as.character(fm$gene_class), var$gene_class)
})

test_that("special characters in obs/var names survive", {
  skip_if_no_anndata()
  n_cells <- 6L; n_genes <- 4L
  X <- matrix(seq_len(n_cells * n_genes) + 0,
              nrow = n_cells, ncol = n_genes)
  rownames(X) <- c("cell-1", "cell.2", "cell_3", "cell:4", "cell 5", "cell/6")
  colnames(X) <- c("ENSG-1", "MT-CO1", "RP/L7", "lncRNA.1")
  ad <- AnnData(X = X)
  s <- convert_anndata_to_seurat(ad)
  expect_equal(colnames(s), rownames(X))
  expect_equal(rownames(s), colnames(X))
})

# ----------------------------------------------------------------------
# C. Layer handling and counts_layer candidate chain
# ----------------------------------------------------------------------

test_that("counts_layer = first match wins when multiple candidates are present", {
  skip_if_no_anndata()
  X <- build_dense_X(10, 6, lambda = 5)
  X2 <- X * 2L
  ad <- AnnData(X = X, layers = list(raw_counts = X, counts = X2))
  s <- convert_anndata_to_seurat(
    ad,
    counts_layer = c("counts", "raw_counts", "raw_count")
  )
  # 'counts' is first in the candidate vector AND present, so it wins.
  expect_equal(
    unname(as.matrix(GetAssayData(s, layer = "counts"))),
    unname(t(X2))
  )
})

test_that("counts_layer falls through to a later candidate", {
  skip_if_no_anndata()
  X <- build_dense_X(10, 6)
  ad <- AnnData(X = X, layers = list(raw_count = X))
  s <- convert_anndata_to_seurat(
    ad,
    counts_layer = c("counts", "raw_counts", "raw_count")
  )
  expect_equal(
    unname(as.matrix(GetAssayData(s, layer = "counts"))),
    unname(t(X))
  )
})

test_that("when no candidate matches, X becomes counts and no separate data layer is added", {
  skip_if_no_anndata()
  X <- build_dense_X(10, 6)
  X_log <- log1p(X)
  ad <- AnnData(X = X, layers = list(logcounts = X_log))
  s <- convert_anndata_to_seurat(ad, counts_layer = c("counts", "raw_counts"))
  # counts == X
  expect_equal(unname(as.matrix(GetAssayData(s, layer = "counts"))),
               unname(t(X)))
  # We did NOT explicitly call SetAssayData(layer = "data") because there
  # was no separate counts layer to fall back to. In Seurat 5 that means the
  # 'data' layer is absent (Layers() omits it).
  expect_false("data" %in% SeuratObject::Layers(s[["RNA"]]))
})

test_that("when a counts layer is present, counts != data and adata$X becomes data", {
  skip_if_no_anndata()
  X <- build_dense_X(10, 6)
  raw <- X * 3L
  ad <- AnnData(X = X, layers = list(counts = raw))
  s <- convert_anndata_to_seurat(ad, counts_layer = "counts")
  counts <- as.matrix(GetAssayData(s, layer = "counts"))
  data <- as.matrix(GetAssayData(s, layer = "data"))
  expect_equal(unname(counts), unname(t(raw)))
  expect_equal(unname(data),   unname(t(X)))
  expect_false(isTRUE(all.equal(unname(counts), unname(data))))
})

test_that("counts_layer = NULL or empty string falls back to X without warning", {
  skip_if_no_anndata()
  X <- build_dense_X(10, 6)
  ad <- AnnData(X = X, layers = list(counts = X * 2L))

  # Capture warnings without suppressing the function's normal status messages.
  captured <- list()
  run <- function(arg) {
    withCallingHandlers(
      suppressMessages(convert_anndata_to_seurat(ad, counts_layer = arg)),
      warning = function(w) {
        captured[[length(captured) + 1L]] <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
  }
  s1 <- run(NULL)
  s2 <- run("")
  expect_length(captured, 0L)

  # Both should use X as counts (since we asked for nothing in particular).
  expect_equal(unname(as.matrix(GetAssayData(s1, layer = "counts"))),
               unname(t(X)))
  expect_equal(unname(as.matrix(GetAssayData(s2, layer = "counts"))),
               unname(t(X)))
})

test_that("Seurat-policy: feature names with underscores get sanitized to dashes (documented)", {
  skip_if_no_anndata()
  # Seurat policy: feature names cannot contain '_'. The conversion still
  # succeeds (the dash-rewritten names match between counts/data/var) and
  # Seurat itself emits the warning. We assert that:
  #   - the warning fires (proving the policy is engaged), and
  #   - the data layer + var feature meta still attach successfully.
  n_cells <- 8L; n_genes <- 4L
  X <- matrix(seq_len(n_cells * n_genes) + 0,
              nrow = n_cells, ncol = n_genes)
  rownames(X) <- sprintf("cell_%02d", seq_len(n_cells))
  colnames(X) <- sprintf("gene_%02d", seq_len(n_genes))
  var <- data.frame(
    gene_class = rep("protein_coding", n_genes),
    row.names = colnames(X),
    stringsAsFactors = FALSE
  )
  ad <- AnnData(X = X, var = var, layers = list(counts = X * 2L))

  warned <- FALSE
  withCallingHandlers(
    s <- suppressMessages(convert_anndata_to_seurat(ad, counts_layer = "counts")),
    warning = function(w) {
      if (grepl("underscore", conditionMessage(w))) warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(warned)
  expect_true(all(grepl("^gene-", rownames(s))))
  # Both counts and data attached cleanly (the bug we fixed).
  expect_true("counts" %in% SeuratObject::Layers(s[["RNA"]]))
  expect_true("data"   %in% SeuratObject::Layers(s[["RNA"]]))
  fm <- s[["RNA"]][[]]
  expect_equal(rownames(fm), rownames(s))
  expect_true("gene_class" %in% colnames(fm))
})

# ----------------------------------------------------------------------
# D. obsm / reductions
# ----------------------------------------------------------------------

test_that("default reduction_map maps standard scanpy keys", {
  skip_if_no_anndata()
  n <- 14
  X <- build_dense_X(n, 8)
  obsm <- list(
    X_pca   = matrix(rnorm(n * 5), n, 5),
    X_umap  = matrix(rnorm(n * 2), n, 2),
    X_tsne  = matrix(rnorm(n * 2), n, 2)
  )
  ad <- AnnData(X = X, obsm = obsm)
  s <- convert_anndata_to_seurat(ad)
  expect_setequal(names(s@reductions), c("pca", "umap", "tsne"))
  expect_equal(s@reductions$pca@key,  "PC_")
  expect_equal(s@reductions$umap@key, "UMAP_")
  expect_equal(s@reductions$tsne@key, "tSNE_")
  expect_equal(unname(s@reductions$pca@cell.embeddings), unname(obsm$X_pca))
})

test_that("non-standard obsm keys get a derived reduction name", {
  skip_if_no_anndata()
  n <- 10
  X <- build_dense_X(n, 5)
  obsm <- list(X_my_emb = matrix(rnorm(n * 4), n, 4))
  ad <- AnnData(X = X, obsm = obsm)
  s <- convert_anndata_to_seurat(ad)
  expect_true("my_emb" %in% names(s@reductions))
  # Embeddings preserved
  expect_equal(unname(s@reductions$my_emb@cell.embeddings), unname(obsm$X_my_emb))
})

test_that("user-supplied reduction_map renames and adds embeddings", {
  skip_if_no_anndata()
  n <- 10
  X <- build_dense_X(n, 5)
  obsm <- list(
    X_pca = matrix(rnorm(n * 3), n, 3),
    X_special = matrix(rnorm(n * 2), n, 2)
  )
  ad <- AnnData(X = X, obsm = obsm)
  s <- convert_anndata_to_seurat(
    ad,
    reduction_map = list(
      X_pca     = list(name = "renamedPCA", key = "rPCA_"),
      X_special = list(name = "special", key = "Sp_")
    )
  )
  expect_true(all(c("renamedPCA", "special") %in% names(s@reductions)))
  expect_equal(s@reductions$renamedPCA@key, "rPCA_")
  expect_equal(s@reductions$special@key, "Sp_")
})

test_that("reduction_map with name=NA disables a default", {
  skip_if_no_anndata()
  n <- 10
  X <- build_dense_X(n, 5)
  obsm <- list(
    X_pca  = matrix(rnorm(n * 3), n, 3),
    X_umap = matrix(rnorm(n * 2), n, 2)
  )
  ad <- AnnData(X = X, obsm = obsm)
  s <- convert_anndata_to_seurat(
    ad,
    reduction_map = list(X_umap = list(name = NA, key = NA))
  )
  expect_true("pca" %in% names(s@reductions))
  expect_false("umap" %in% names(s@reductions))
})

# ----------------------------------------------------------------------
# E. uns / orig.ident logic
# ----------------------------------------------------------------------

test_that("uns$conversion_source becomes orig.ident when no obs$orig.ident", {
  skip_if_no_anndata()
  X <- build_dense_X(8, 5)
  ad <- AnnData(X = X, uns = list(conversion_source = "fixture42"))
  s <- convert_anndata_to_seurat(ad)
  expect_equal(unique(as.character(s@meta.data$orig.ident)), "fixture42")
})

test_that("uns with arbitrary unrelated keys does not break", {
  skip_if_no_anndata()
  X <- build_dense_X(8, 5)
  ad <- AnnData(X = X, uns = list(some_param = 42, notes = "hello"))
  s <- convert_anndata_to_seurat(ad)
  expect_equal(unique(as.character(s@meta.data$orig.ident)), "AnnData")
})

test_that("explicit orig.ident parameter beats uns and file basename", {
  skip_if_no_anndata()
  X <- build_dense_X(6, 4)
  ad <- AnnData(X = X, uns = list(conversion_source = "from_uns"))
  s <- convert_anndata_to_seurat(ad, orig.ident = "explicit_value")
  expect_equal(unique(as.character(s@meta.data$orig.ident)), "explicit_value")
})

# ----------------------------------------------------------------------
# F. File-path entry
# ----------------------------------------------------------------------

test_that("path entry: orig.ident defaults to file basename when no uns", {
  skip_if_no_anndata()
  X <- build_dense_X(8, 4)
  ad <- AnnData(X = X)
  path <- file.path(tempdir(), "my_sample_42.h5ad")
  on.exit(unlink(path), add = TRUE)
  write_h5ad(ad, path)
  s <- convert_anndata_to_seurat(path)
  expect_equal(unique(as.character(s@meta.data$orig.ident)), "my_sample_42")
})

test_that("path entry: tilde expansion works", {
  skip_if_no_anndata()
  X <- build_dense_X(6, 4)
  ad <- AnnData(X = X)
  rel <- file.path("~", "convert2anndata_test_tilde.h5ad")
  abs <- path.expand(rel)
  on.exit(unlink(abs), add = TRUE)
  write_h5ad(ad, abs)
  s <- convert_anndata_to_seurat(rel)
  expect_s4_class(s, "Seurat")
})

test_that("path entry: nonexistent file gives a clear error", {
  expect_error(
    convert_anndata_to_seurat("/nope/does/not/exist.h5ad"),
    regexp = "does not exist"
  )
})

test_that("path entry: invalid h5ad file gives a wrapped error", {
  skip_if_no_anndata()
  bogus <- tempfile(fileext = ".h5ad")
  writeLines("this is not an h5ad file", bogus)
  on.exit(unlink(bogus), add = TRUE)
  expect_error(
    convert_anndata_to_seurat(bogus),
    regexp = "Failed to read AnnData h5ad file"
  )
})

# ----------------------------------------------------------------------
# G. Custom assay name
# ----------------------------------------------------------------------

test_that("assay name can be overridden and the data layer attaches to it", {
  skip_if_no_anndata()
  X <- build_dense_X(10, 6)
  raw <- X * 2L
  ad <- AnnData(X = X, layers = list(counts = raw))
  s <- convert_anndata_to_seurat(ad, counts_layer = "counts", assay = "ATAC")
  # Use names(s@assays) instead of Seurat::Assays() because the latter
  # generic gets shadowed by SummarizedExperiment in cross-file test_dir runs.
  expect_true("ATAC" %in% names(s@assays))
  expect_false("RNA" %in% names(s@assays))
  counts <- as.matrix(GetAssayData(s, assay = "ATAC", layer = "counts"))
  data   <- as.matrix(GetAssayData(s, assay = "ATAC", layer = "data"))
  expect_equal(unname(counts), unname(t(raw)))
  expect_equal(unname(data),   unname(t(X)))
})

test_that("reductions get the custom assay name when assay is overridden", {
  skip_if_no_anndata()
  n <- 10
  X <- build_dense_X(n, 6)
  ad <- AnnData(X = X, obsm = list(X_pca = matrix(rnorm(n * 3), n, 3)))
  s <- convert_anndata_to_seurat(ad, assay = "ATAC")
  expect_equal(s@reductions$pca@assay.used, "ATAC")
})

# ----------------------------------------------------------------------
# H. Round-trip integrity
# ----------------------------------------------------------------------

test_that("AnnData -> Seurat -> SCE -> AnnData -> Seurat preserves counts and obs", {
  skip_if_no_anndata()
  n <- 16; m <- 10
  X <- build_dense_X(n, m, lambda = 2)
  obs <- data.frame(
    cell_type = factor(rep(c("A", "B"), n / 2)),
    n_counts  = rowSums(X),
    row.names = rownames(X)
  )
  ad <- AnnData(X = X, obs = obs,
                layers = list(counts = X),
                obsm = list(X_pca = matrix(rnorm(n * 4), n, 4)))

  s1 <- convert_anndata_to_seurat(ad, counts_layer = "counts")
  sce <- convert_seurat_to_sce(s1)
  ad2 <- convert_to_anndata(sce, assayName = "counts")
  s2 <- convert_anndata_to_seurat(ad2, counts_layer = "counts")

  # Counts preserved
  expect_equal(
    unname(as.matrix(GetAssayData(s1, layer = "counts"))),
    unname(as.matrix(GetAssayData(s2, layer = "counts")))
  )
  # cell_type column survives
  expect_equal(
    as.character(s1@meta.data$cell_type),
    as.character(s2@meta.data$cell_type)
  )
  # obs/var names preserved
  expect_equal(colnames(s1), colnames(s2))
  expect_equal(rownames(s1), rownames(s2))
})

# ----------------------------------------------------------------------
# I. Robustness
# ----------------------------------------------------------------------

test_that("empty obsm produces zero reductions", {
  skip_if_no_anndata()
  X <- build_dense_X(8, 5)
  ad <- AnnData(X = X)
  s <- convert_anndata_to_seurat(ad)
  expect_length(s@reductions, 0)
})

test_that("obsm with mismatched cell count is skipped (warn) but valid ones attach", {
  skip_if_no_anndata()
  n <- 8
  X <- build_dense_X(n, 5)
  good <- matrix(rnorm(n * 2), n, 2)
  bad  <- matrix(rnorm((n - 1) * 2), n - 1, 2)
  ad <- AnnData(X = X, obsm = list(X_pca = good))

  # extract_anndata_obsm + attach_reductions_seurat is the public path; bad
  # obsm is normally rejected by AnnData itself, so we feed it directly to
  # attach_reductions_seurat to verify the guard.
  s <- convert_anndata_to_seurat(ad)
  bad_named <- list(X_pca = good, X_bad = bad)
  expect_warning(
    s2 <- attach_reductions_seurat(s, bad_named),
    regexp = "skipping"
  )
  expect_true("pca" %in% names(s2@reductions))
})

test_that("very small AnnData (3x2) and skinny (2x100) do not break the converter", {
  skip_if_no_anndata()
  ad_small <- AnnData(X = matrix(c(1, 2, 3, 4, 5, 6), 3, 2))
  s1 <- convert_anndata_to_seurat(ad_small)
  expect_equal(ncol(s1), 3); expect_equal(nrow(s1), 2)

  ad_skinny <- AnnData(X = matrix(rpois(200, 1), 2, 100))
  s2 <- convert_anndata_to_seurat(ad_skinny)
  expect_equal(ncol(s2), 2); expect_equal(nrow(s2), 100)
})

test_that("standalone shim and package function produce identical Seurat objects", {
  skip_if_no_anndata()
  src <- "/fh/fast/setty_m/user/yhuang2/sarah-nexus/work/convert2seurat/convert_anndata2seurat.R"
  skip_if_not(file.exists(src), "standalone shim not available")
  source(src, local = TRUE)

  X <- build_dense_X(20, 8)
  ad <- AnnData(X = X, layers = list(counts = X * 2L),
                obsm = list(X_pca = matrix(rnorm(60), 20, 3)))
  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  write_h5ad(ad, path)

  s_shim <- convert_anndata2seurat(path, raw_layer_name = "counts")
  s_pkg <- convert_anndata_to_seurat(path, counts_layer = "counts")

  expect_equal(dim(s_shim), dim(s_pkg))
  expect_equal(
    unname(as.matrix(GetAssayData(s_shim, layer = "counts"))),
    unname(as.matrix(GetAssayData(s_pkg, layer = "counts")))
  )
  expect_equal(
    unname(as.matrix(GetAssayData(s_shim, layer = "data"))),
    unname(as.matrix(GetAssayData(s_pkg, layer = "data")))
  )
  expect_equal(sort(names(s_shim@reductions)), sort(names(s_pkg@reductions)))
})

# ----------------------------------------------------------------------
# J. orig.ident precedence laddered
# ----------------------------------------------------------------------

test_that("orig.ident precedence: obs column > arg > uns > path-basename > 'AnnData'", {
  skip_if_no_anndata()
  X <- build_dense_X(6, 4)

  # 1. obs.orig.ident wins over everything (we don't override an existing column).
  obs1 <- data.frame(orig.ident = paste0("from_obs_", 1:6),
                     row.names = rownames(X))
  ad1 <- AnnData(X = X, obs = obs1, uns = list(conversion_source = "from_uns"))
  path1 <- tempfile(fileext = ".h5ad"); on.exit(unlink(path1), add = TRUE)
  write_h5ad(ad1, path1)
  s1 <- convert_anndata_to_seurat(path1, orig.ident = "from_arg")
  expect_equal(s1@meta.data$orig.ident, paste0("from_obs_", 1:6))

  # 2. With no obs.orig.ident, the explicit arg wins over uns and path.
  ad2 <- AnnData(X = X, uns = list(conversion_source = "from_uns"))
  path2 <- tempfile(fileext = ".h5ad"); on.exit(unlink(path2), add = TRUE)
  write_h5ad(ad2, path2)
  s2 <- convert_anndata_to_seurat(path2, orig.ident = "from_arg")
  expect_equal(unique(s2@meta.data$orig.ident), "from_arg")

  # 3. With no obs and no arg, uns wins over path.
  s3 <- convert_anndata_to_seurat(path2)
  expect_equal(unique(s3@meta.data$orig.ident), "from_uns")

  # 4. Without uns either, path basename wins.
  ad4 <- AnnData(X = X)
  path4 <- file.path(tempdir(), "named_sample.h5ad"); on.exit(unlink(path4), add = TRUE)
  write_h5ad(ad4, path4)
  s4 <- convert_anndata_to_seurat(path4)
  expect_equal(unique(s4@meta.data$orig.ident), "named_sample")

  # 5. No obs, no uns, no path -> "AnnData".
  s5 <- convert_anndata_to_seurat(ad4)
  expect_equal(unique(s5@meta.data$orig.ident), "AnnData")
})
