library(testthat)
library(convert2anndata)
library(Matrix)
library(anndata)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)

# ======================================================================
# Complex Seurat object round-trip fidelity (issue #18).
#
#   Seurat -> SCE -> AnnData -> Seurat
#
# (there is no direct Seurat<->AnnData; the cycle goes through SCE both ways).
# COMPLEX Seurat: processed RNA assay (normalised/scaled/PCA), a UMAP
# reduction, nearest-neighbour graphs, a second ADT assay, and factor + NA
# cell metadata. Fixture in helper-roundtrip.R.
# ======================================================================

TOL_F32 <- 1e-4

# Full Seurat -> SCE -> AnnData -> Seurat cycle (in memory).
rt_seu <- function(seu) {
  sce <- rt_quiet(convert_seurat_to_sce(seu))
  ad  <- rt_quiet(convert_to_anndata(sce, assayName = "RNA_counts", useAltExp = TRUE))
  rt_quiet(convert_anndata_to_seurat(ad, counts_layer = "RNA_counts"))
}

# ----------------------------------------------------------------------
# Seurat -> SCE -> AnnData -> Seurat
# ----------------------------------------------------------------------

test_that("complex Seurat round-trips RNA counts and dimensions", {
  skip_if_no_anndata()
  seu  <- build_complex_seurat()
  back <- rt_seu(seu)

  expect_s4_class(back, "Seurat")
  expect_equal(ncol(back), ncol(seu))   # cells
  expect_equal(nrow(back), nrow(seu))   # RNA features
  c0 <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts"))
  c1 <- as.matrix(GetAssayData(back, layer = "counts"))
  expect_lte(rt_max_abs_diff(c0, c1), TOL_F32)
})

test_that("complex Seurat round-trips reductions (pca + umap)", {
  skip_if_no_anndata()
  seu  <- build_complex_seurat()
  back <- rt_seu(seu)

  expect_true(all(c("pca", "umap") %in% names(back@reductions)))
  expect_lte(rt_max_abs_diff(Embeddings(seu, "pca"),
                             Embeddings(back, "pca")), TOL_F32)
  expect_lte(rt_max_abs_diff(Embeddings(seu, "umap"),
                             Embeddings(back, "umap")), TOL_F32)
})

test_that("complex Seurat round-trips cell metadata (factor + NA)", {
  skip_if_no_anndata()
  seu  <- build_complex_seurat()
  back <- rt_seu(seu)

  expect_true(all(c("cell_type", "qc_na") %in% colnames(back@meta.data)))
  expect_equal(as.character(back@meta.data$cell_type),
               as.character(seu@meta.data$cell_type))
  expect_true(is.na(back@meta.data$qc_na[1]))
})

test_that("complex Seurat round-trips NN graphs as connectivity", {
  skip_if_no_anndata()
  # Seurat Graph -> SCE colPair -> AnnData obsp -> Seurat Graph. Edge weights
  # are dropped by convert_graph_to_colPair (binarised), and the obsp key gains
  # the assay prefix (RNA_nn -> RNA_RNA_nn), so we assert that *some* returned
  # graph reproduces the source nn connectivity pattern.
  seu  <- build_complex_seurat()
  back <- rt_seu(seu)

  expect_gt(length(back@graphs), 0)
  src <- rt_pattern(seu@graphs[["RNA_nn"]])
  matches <- vapply(back@graphs, function(g) {
    p <- rt_pattern(g)
    is.matrix(p) && all(dim(p) == dim(src)) && all(p == src)
  }, logical(1))
  expect_true(any(matches),
              info = "no round-tripped graph matched the source RNA_nn pattern")
})

# ----------------------------------------------------------------------
# Seurat -> SCE: a secondary assay becomes an altExp
# ----------------------------------------------------------------------

test_that("a secondary Seurat assay (ADT) becomes an SCE altExp", {
  skip_if_no_anndata()
  seu <- build_complex_seurat()
  sce <- rt_quiet(convert_seurat_to_sce(seu))
  expect_true("ADT" %in% altExpNames(sce))
  expect_equal(nrow(altExp(sce, "ADT")), 6L)  # 6 ADT features
})

# ----------------------------------------------------------------------
# Known gap (documented, not silently missing) -- see TODO(#18)
# ----------------------------------------------------------------------

test_that("[known gap] a secondary Seurat assay is not restored as a Seurat assay", {
  # TODO(#18): a 2nd Seurat assay -> SCE altExp -> AnnData uns$altExperiments,
  # but convert_anndata_to_seurat() does not rebuild Seurat assays from
  # uns$altExperiments, so ADT does not return as a Seurat assay across the
  # full cycle.
  skip("Known gap: convert_anndata_to_seurat() does not restore altExperiments.")
  seu  <- build_complex_seurat()
  back <- rt_seu(seu)
  expect_true("ADT" %in% names(back@assays))
})
