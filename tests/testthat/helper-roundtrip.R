# Shared scaffolding for the complex round-trip test files
# (test-roundtrip_complex_{anndata,sce,seurat}.R).
#
# Reuses the same skip guard and seeded-RNG idiom already established across
# the suite (test-comprehensive_anndata_to_seurat.R, test-raw_and_obsp.R,
# test-realistic_scanpy_dataset.R), centralised here so the three round-trip
# files do not each rebuild the (large) complex fixtures. testthat sources
# helper-*.R before the test files, so these are visible everywhere.

# ---- Guards / RNG (mirror the existing in-file definitions) ----------------

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
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv))
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)))
  }
  set.seed(seed)
  force(expr)
}

# ---- Small utilities -------------------------------------------------------

# The conversion functions emit a lot of timestamped status via cat(); swallow
# stdout but let warnings/errors surface so the test log stays readable.
rt_quiet <- function(expr) {
  utils::capture.output(res <- force(expr))
  res
}

# Robust AnnData mapping keys (layers/obsm/obsp/varm/varp/uns) -> character.
rt_keys <- function(mapping) {
  bi <- reticulate::import_builtins()
  tryCatch(as.character(bi$list(mapping$keys())),
           error = function(e) as.character(names(mapping)))
}

# scanpy obsm keys carry an "X_" prefix; reducedDimNames do not. Normalise so
# X_pca / PCA / pca all compare equal.
rt_norm_key <- function(k) tolower(sub("^X_", "", k))

rt_max_abs_diff <- function(a, b) {
  a <- as.matrix(a); b <- as.matrix(b)
  if (!all(dim(a) == dim(b))) return(Inf)
  max(abs(unname(a) - unname(b)))
}

# Logical sparsity pattern (which entries are non-zero), for graph/obsp
# comparisons where edge *weights* may be dropped but connectivity must hold.
rt_pattern <- function(m) {
  m <- as.matrix(m)
  unname(m != 0)
}

# Full AnnData round trip through an on-disk .h5ad (the realistic path).
rt_write_read <- function(adata_obj) {
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  anndata::write_h5ad(adata_obj, tmp)
  anndata::read_h5ad(tmp)
}

# ---- Complex fixtures ------------------------------------------------------

# A COMPLEX AnnData exercising every component the converters touch:
# sparse OR dense X, two layers (counts/logcounts), two obsm reductions,
# varm gene loadings, two obsp graphs, a varp graph, raw, and obs/var with
# factor / NA / logical / character / numeric columns.
build_complex_anndata <- function(n_cells = 24L, n_genes = 12L,
                                  sparse_X = TRUE, with_raw = TRUE, seed = 1L) {
  with_seed(seed, {
    np <- reticulate::import("numpy")
    sp <- reticulate::import("scipy.sparse")

    X <- matrix(rpois(n_cells * n_genes, lambda = 3), n_cells, n_genes)
    rownames(X) <- sprintf("cell-%02d", seq_len(n_cells))
    colnames(X) <- sprintf("gene-%02d", seq_len(n_genes))
    logX <- log1p(X)

    obs <- data.frame(
      cell_type = factor(rep(c("A", "B", "C"), length.out = n_cells),
                         levels = c("A", "B", "C")),
      n_counts  = rowSums(X),
      pct_mt    = c(NA_real_, runif(n_cells - 1L, 0, 10)),
      is_ctrl   = rep(c(TRUE, FALSE), length.out = n_cells),
      donor     = sprintf("D%d", rep(1:2, length.out = n_cells)),
      row.names = rownames(X), stringsAsFactors = FALSE
    )
    var <- data.frame(
      gene_class      = factor(rep(c("pc", "lnc"), length.out = n_genes)),
      highly_variable = rep(c(TRUE, FALSE), length.out = n_genes),
      mean_expr       = colMeans(X),
      row.names = colnames(X), stringsAsFactors = FALSE
    )

    conn <- abs(Matrix::rsparsematrix(n_cells, n_cells,
                                      density = 0.3, symmetric = TRUE))
    dimnames(conn) <- list(rownames(X), rownames(X))

    X_in <- if (sparse_X) sp$csr_matrix(np$asarray(logX, dtype = "float32")) else logX

    ad <- anndata::AnnData(
      X      = X_in,
      obs    = obs,
      var    = var,
      layers = list(counts = X, logcounts = logX),
      obsm   = list(X_pca  = matrix(rnorm(n_cells * 5), n_cells, 5),
                    X_umap = matrix(rnorm(n_cells * 2), n_cells, 2)),
      varm   = list(PCs = matrix(rnorm(n_genes * 5), n_genes, 5)),
      obsp   = list(connectivities = conn, distances = conn * 2),
      varp   = list(corr = abs(Matrix::rsparsematrix(n_genes, n_genes,
                                                     density = 0.4, symmetric = TRUE))),
      uns    = list(conversion_source = "complex_fixture")
    )
    if (with_raw) ad$raw <- anndata::AnnData(X = X)  # same gene set raw counts
    ad
  })
}

# A COMPLEX SingleCellExperiment: three assays (sparse counts + dense
# logcounts + dense scaledata), two reducedDims, an altExp (multimodal ADT),
# a colPair kNN graph, colData (factor/NA/mixed) and rowData metadata.
build_complex_sce <- function(n_cells = 30L, n_genes = 15L, seed = 2L) {
  with_seed(seed, {
    cnt <- abs(Matrix::rsparsematrix(n_genes, n_cells, density = 0.5))
    dimnames(cnt) <- list(sprintf("g%02d", seq_len(n_genes)),
                          sprintf("c%02d", seq_len(n_cells)))
    dense  <- as.matrix(cnt)
    logc   <- log1p(dense)                                 # dense
    # A third, distinct dense assay. Deterministic and NaN-free (avoids the
    # zero-variance-gene trap that scale() would hit on a sparse count matrix).
    scaled <- (dense - rowMeans(dense)) / (apply(dense, 1, sd) + 1)

    cd <- S4Vectors::DataFrame(
      cell_type = factor(rep(c("X", "Y"), length.out = n_cells)),
      qc        = c(NA_real_, runif(n_cells - 1L)),
      batch     = rep(c("b1", "b2", "b3"), length.out = n_cells),
      row.names = colnames(cnt)
    )
    rd <- S4Vectors::DataFrame(
      symbol = sprintf("SYM%d", seq_len(n_genes)),
      is_hvg = rep(c(TRUE, FALSE), length.out = n_genes),
      row.names = rownames(cnt)
    )
    pca  <- matrix(rnorm(n_cells * 6), n_cells, 6, dimnames = list(colnames(cnt), NULL))
    umap <- matrix(rnorm(n_cells * 2), n_cells, 2, dimnames = list(colnames(cnt), NULL))

    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = cnt, logcounts = logc, scaledata = scaled),
      colData = cd, rowData = rd,
      reducedDims = list(PCA = pca, UMAP = umap)
    )

    adt <- matrix(abs(rnorm(4 * n_cells)), 4, n_cells,
                  dimnames = list(sprintf("ADT%d", 1:4), colnames(cnt)))
    SingleCellExperiment::altExp(sce, "ADT") <-
      SingleCellExperiment::SingleCellExperiment(assays = list(counts = adt))

    knn <- BiocNeighbors::findKNN(pca, k = 5)
    SingleCellExperiment::colPair(sce, "knn") <- S4Vectors::SelfHits(
      from = rep(seq_len(n_cells), each = 5),
      to   = as.vector(t(knn$index)),
      x    = as.vector(t(knn$distance)),
      nnode = n_cells
    )
    S4Vectors::metadata(sce)$run_info <- "complex-sce"
    sce
  })
}

# A COMPLEX Seurat object: a processed RNA assay (normalised + scaled + PCA),
# a manually-attached UMAP reduction (avoids the uwot dependency / nondeterminism),
# nearest-neighbour graphs, a second ADT assay, and factor + NA cell metadata.
build_complex_seurat <- function(n_genes = 40L, n_cells = 60L, seed = 3L) {
  with_seed(seed, {
    cm <- matrix(rpois(n_genes * n_cells, lambda = 3), n_genes, n_cells,
                 dimnames = list(sprintf("g%02d", seq_len(n_genes)),
                                 sprintf("c%02d", seq_len(n_cells))))
    seu <- Seurat::CreateSeuratObject(counts = cm)
    seu <- Seurat::NormalizeData(seu, verbose = FALSE)
    seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
    seu <- Seurat::ScaleData(seu, verbose = FALSE)
    seu <- Seurat::RunPCA(seu, npcs = 10, verbose = FALSE)
    seu <- Seurat::FindNeighbors(seu, dims = 1:10, verbose = FALSE)

    # Manual UMAP reduction from the PCA embeddings: deterministic, no uwot.
    emb <- Seurat::Embeddings(seu, "pca")[, 1:2, drop = FALSE]
    colnames(emb) <- c("UMAP_1", "UMAP_2")
    seu[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = emb, key = "UMAP_", assay = "RNA"
    )

    seu$cell_type <- factor(rep(c("T", "B", "NK"), length.out = n_cells))
    seu$qc_na <- c(NA_real_, runif(n_cells - 1L))

    adt <- matrix(rpois(6 * n_cells, lambda = 5), 6, n_cells,
                  dimnames = list(sprintf("adt%d", 1:6), colnames(seu)))
    seu[["ADT"]] <- Seurat::CreateAssayObject(counts = adt)
    seu
  })
}
