library(testthat)
library(convert2anndata)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(S4Vectors)
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

mock_h5ad_path <- function() {
  skip_if_no_anndata()
  n_cells <- 20
  n_genes <- 30
  X <- matrix(rpois(n_cells * n_genes, lambda = 3), nrow = n_cells, ncol = n_genes)
  rownames(X) <- paste0("cell", seq_len(n_cells))
  colnames(X) <- paste0("gene", seq_len(n_genes))

  obs <- data.frame(
    cell_type = rep(c("A", "B"), n_cells / 2),
    n_counts = rowSums(X),
    row.names = rownames(X)
  )
  var <- data.frame(
    gene_type = rep(c("g1", "g2"), n_genes / 2),
    row.names = colnames(X)
  )
  pca <- matrix(runif(n_cells * 3), nrow = n_cells, ncol = 3)
  rownames(pca) <- rownames(X)

  adata <- AnnData(
    X = X,
    obs = obs,
    var = var,
    layers = list(counts = X),
    obsm = list(X_pca = pca, X_umap = pca[, 1:2])
  )
  path <- tempfile(fileext = ".h5ad")
  write_h5ad(adata, path)
  path
}

test_that("convert_anndata_to_seurat preserves dims, metadata, and reductions", {
  skip_if_not_installed("Seurat")
  path <- mock_h5ad_path()
  adata <- read_h5ad(path)

  seurat_obj <- convert_anndata_to_seurat(adata, counts_layer = "counts")

  expect_s4_class(seurat_obj, "Seurat")
  expect_equal(ncol(seurat_obj), nrow(adata))
  expect_equal(nrow(seurat_obj), ncol(adata))
  expect_equal(colnames(seurat_obj), as.character(adata$obs_names))
  expect_equal(rownames(seurat_obj), as.character(adata$var_names))

  expect_true("cell_type" %in% colnames(seurat_obj@meta.data))
  expect_equal(
    as.character(seurat_obj@meta.data$cell_type),
    as.character(adata$obs$cell_type)
  )

  expect_true("pca" %in% names(seurat_obj@reductions))
  expect_true("umap" %in% names(seurat_obj@reductions))
  expect_equal(
    unname(seurat_obj@reductions$pca@cell.embeddings),
    unname(as.matrix(adata$obsm[["X_pca"]]))
  )
})

test_that("convert_anndata_to_sce preserves dims and reducedDims", {
  path <- mock_h5ad_path()
  adata <- read_h5ad(path)

  sce <- convert_anndata_to_sce(adata, counts_layer = "counts")

  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(ncol(sce), nrow(adata))
  expect_equal(nrow(sce), ncol(adata))
  expect_true("PCA" %in% names(reducedDims(sce)))
  expect_equal(
    unname(reducedDim(sce, "PCA")),
    unname(as.matrix(adata$obsm[["X_pca"]]))
  )
  expect_equal(as.character(colData(sce)$cell_type), as.character(adata$obs$cell_type))
})

test_that("anndata -> seurat -> anndata roundtrip preserves counts", {
  skip_if_not_installed("Seurat")
  path <- mock_h5ad_path()
  adata <- read_h5ad(path)

  seurat_obj <- convert_anndata_to_seurat(adata, counts_layer = "counts")
  sce_back <- convert_seurat_to_sce(seurat_obj)
  adata_back <- convert_to_anndata(sce_back, assayName = "counts")

  expect_equal(dim(adata_back$X), dim(adata$X))
  expect_equal(
    unname(as.matrix(adata_back$X)),
    unname(as.matrix(adata$X))
  )
})

test_that("convert_anndata_to_seurat handles missing counts layer", {
  skip_if_not_installed("Seurat")
  skip_if_no_anndata()
  n_cells <- 10
  n_genes <- 12
  X <- matrix(rpois(n_cells * n_genes, lambda = 2), nrow = n_cells, ncol = n_genes)
  rownames(X) <- paste0("c", seq_len(n_cells))
  colnames(X) <- paste0("g", seq_len(n_genes))
  adata <- AnnData(X = X)

  seurat_obj <- convert_anndata_to_seurat(adata, counts_layer = "counts")

  expect_s4_class(seurat_obj, "Seurat")
  expect_equal(ncol(seurat_obj), n_cells)
  expect_equal(nrow(seurat_obj), n_genes)
})
