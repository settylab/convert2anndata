library(testthat)
library(convert2anndata)
library(SingleCellExperiment)
library(S4Vectors)
library(Matrix)
library(SummarizedExperiment)
library(anndata)
library(R6)
library(BiocNeighbors)


create_mock_sce <- function() {
  n_cells <- 50
  n_genes <- 100

  # Create a main assay matrix
  counts <- matrix(runif(n_cells * n_genes), nrow = n_genes, ncol = n_cells)

  # Create an alternative experiment with matching column numbers
  alt_counts <- matrix(runif(5 * n_cells), nrow = 5, ncol = n_cells)
  alt1 <- SingleCellExperiment(assays = list(counts = alt_counts))

  # Create reduced dimensions
  pca <- matrix(runif(3 * n_cells), nrow = n_cells, ncol = 3)
  reducedDims <- list(PCA = pca, UMAP = pca, other_rep = pca)

  # Create colData and rowData
  colData <- DataFrame(cell_type = rep(c("A", "B"), n_cells / 2))
  rowData <- DataFrame(gene_type = rep(c("gene1", "gene2"), n_genes / 2))

  # Create the main SCE object
  sce <- SingleCellExperiment(
    assays = list(counts = counts, logcounts = counts),
    colData = colData,
    rowData = rowData,
    reducedDims = reducedDims,
    altExps = list(alt1 = alt1)
  )

  # Compute k-nearest neighbors and add to the SCE object
  knn <- findKNN(pca, k = 10)
  colPairs(sce) <- knn

  # Add additional metadata
  metadata(sce) <- list(experiment_info = "Mock experiment")

  return(sce)
}

test_that("convert_to_anndata works with Seurat object and preserves commands in a serializable format", {
  skip_if_not_installed("Seurat")
  library(Seurat)

  # Create Seurat object with command history
  seurat_obj <- CreateSeuratObject(counts = matrix(rpois(200, lambda = 5), nrow = 20, ncol = 10))
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 5)

  # Convert to SCE
  sce <- convert_seurat_to_sce(seurat_obj)

  # Convert to AnnData
  ad <- convert_to_anndata(sce)

  # Ensure it is an R6 object
  expect_true(R6::is.R6(ad))

  # Check that dimensions are preserved
  expect_equal(dim(ad$X), rev(dim(sce)))

  # Ensure commands are present in `uns`
  expect_true("commands" %in% names(ad$uns))
  expect_true(!is.null(ad$uns$commands))

  # Ensure commands are **fully serializable**
  expect_true(all(sapply(ad$uns$commands, function(cmd) {
    all(sapply(cmd, function(value) {
      is.null(value) || is.atomic(value) || is.list(value)
    }))
  })))

  # Ensure that writing to HDF5 does not raise errors
  temp_h5ad <- tempfile(fileext = ".h5ad")
  expect_silent(write_h5ad(ad, temp_h5ad))

  # Ensure written H5AD file is readable
  expect_silent(read_h5ad(temp_h5ad))
})

test_that("convert_to_anndata works correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  tryCatch(
    {
      ad <- convert_to_anndata(sce)
      expect_true(R6::is.R6(ad)) # Ensure it is an R6 object
      expect_equal(nrow(ad$X), 2)
      expect_equal(ncol(ad$X), 2)
      if (!is.null(ad$obs)) {
        expect_equal(nrow(ad$obs), 2)
      }
      if (!is.null(ad$var)) {
        expect_equal(nrow(ad$var), 2)
      }
    },
    error = function(e) {
      print(reticulate::py_last_error())
      stop(e)
    }
  )
})

test_that("convert_to_anndata handles altExps correctly", {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  alt_sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(5:8, ncol = 2)))
  altExp(sce, "alt1") <- alt_sce
  tryCatch(
    {
      ad <- convert_to_anndata(sce, useAltExp = TRUE)
      expect_true("altExperiments" %in% names(ad$uns))
      expect_true(R6::is.R6(ad$uns$altExperiments$alt1)) # Ensure it is an R6 object
      if (!is.null(ad$obs)) {
        expect_equal(nrow(ad$obs), 2)
      }
      if (!is.null(ad$var)) {
        expect_equal(nrow(ad$var), 2)
      }
    },
    error = function(e) {
      print(reticulate::py_last_error())
      stop(e)
    }
  )
})

test_that("convert_to_anndata works with Seurat objects", {
  skip_if_not_installed("Seurat")
  library(Seurat)
  seurat_obj <- Seurat::CreateSeuratObject(counts = matrix(1:4, ncol = 2))
  sce <- convert_seurat_to_sce(seurat_obj)
  ad <- convert_to_anndata(sce)
  expect_true(R6::is.R6(ad)) # Ensure it is an R6 object
  expect_equal(nrow(ad$X), 2)
  expect_equal(ncol(ad$X), 2)
  if (!is.null(ad$obs)) {
    expect_equal(nrow(ad$obs), 2)
  }
  if (!is.null(ad$var)) {
    expect_equal(nrow(ad$var), 2)
  }
})

test_that("convert_to_anndata works with complex SCE input", {
  sce <- create_mock_sce()

  # Run the conversion function
  ad <- convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE)

  # Check the main data matrix
  mat1 <- as.matrix(ad$X)
  mat2 <- t(as.matrix(assay(sce, "counts")))
  # Remove dimnames attributes before comparison
  mat1 <- unname(mat1)
  mat2 <- unname(mat2)
  # Compare the matrices
  expect_equal(dim(mat1), dim(mat2))
  expect_equal(mat1, mat2, tolerance = 1e-6)

  # Check the layers (if you have multiple assays)
  expect_true("logcounts" %in% names(ad$layers))
  mat1 <- as.matrix(ad$layers[["logcounts"]])
  mat2 <- t(as.matrix(assay(sce, "logcounts")))
  # Remove dimnames attributes before comparison
  mat1 <- unname(mat1)
  mat2 <- unname(mat2)
  # Compare the matrices
  expect_equal(dim(mat1), dim(mat2))
  expect_equal(mat1, mat2, tolerance = 1e-6)

  # Check the dimensional reductions
  expect_true("X_pca" %in% names(ad$obsm))
  expect_equal(ad$obsm$X_pca, reducedDim(sce, "PCA"))

  # Check the obs data
  expect_equal(ad$obs$cell_type, colData(sce)$cell_type)

  # Check the var data
  expect_equal(ad$var$gene_type, rowData(sce)$gene_type)

  # Check the altExp (if present)
  alt_exp_name <- "altExperiments"
  expect_true(alt_exp_name %in% names(ad$uns))
  expect_true("alt1" %in% names(ad$uns[[alt_exp_name]]))

  mat1 <- as.matrix(ad$uns[[alt_exp_name]]$alt1$X)
  mat2 <- t(as.matrix(assay(altExp(sce, "alt1"), "counts")))
  # Remove dimnames attributes before comparison
  mat1 <- unname(mat1)
  mat2 <- unname(mat2)
  # Compare the matrices
  expect_equal(mat1, mat2, tolerance = 1e-6)

  # Check the colPairs data (if present)
  col_pairs <- c("distance", "connectivity")
  for (name in names(col_pairs)) {
    expect_true(name %in% names(ad$obsp))
  }
})
