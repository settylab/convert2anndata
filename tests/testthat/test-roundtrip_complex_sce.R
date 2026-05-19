library(testthat)
library(convert2anndata)
library(Matrix)
library(anndata)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(BiocNeighbors)

# ======================================================================
# Complex SingleCellExperiment round-trip fidelity (issue #18).
#
#   SCE -> AnnData -> SCE
#
# on a COMPLEX SCE: sparse counts + two dense assays, PCA + UMAP reducedDims,
# an altExp (multimodal ADT), a colPair kNN graph, and colData/rowData with
# factor / NA / character columns. Fixture in helper-roundtrip.R.
#
# The round trip is done in memory (not through .h5ad): convert_to_anndata()
# stores altExps as nested AnnData objects under uns$altExperiments, which is
# not a goal of HDF5 serialisation. Components carried only one way are
# documented as skip()'d tests with TODO(#18).
# ======================================================================

TOL_F32 <- 1e-4   # counts pass through float32 inside AnnData

rt_sce <- function(sce, assayName = "counts") {
  ad <- rt_quiet(convert_to_anndata(sce, assayName = assayName, useAltExp = TRUE))
  rt_quiet(convert_anndata_to_sce(ad, counts_layer = "counts"))
}

# ----------------------------------------------------------------------
# SCE -> AnnData -> SCE
# ----------------------------------------------------------------------

test_that("complex SCE round-trips all assays (sparse counts + dense layers)", {
  skip_if_no_anndata()
  skip_if_not_installed("BiocNeighbors")
  sce  <- build_complex_sce()
  back <- rt_sce(sce)

  expect_true(all(c("counts", "logcounts", "scaledata") %in% assayNames(back)))
  expect_lte(rt_max_abs_diff(assay(sce, "counts"), assay(back, "counts")), TOL_F32)
  expect_lte(rt_max_abs_diff(assay(sce, "logcounts"), assay(back, "logcounts")), TOL_F32)
  expect_lte(rt_max_abs_diff(assay(sce, "scaledata"), assay(back, "scaledata")), TOL_F32)
})

test_that("complex SCE round-trips reducedDims (PCA + UMAP)", {
  skip_if_no_anndata()
  skip_if_not_installed("BiocNeighbors")
  sce  <- build_complex_sce()
  back <- rt_sce(sce)

  expect_true(all(c("PCA", "UMAP") %in% reducedDimNames(back)))
  expect_lte(rt_max_abs_diff(reducedDim(sce, "PCA"), reducedDim(back, "PCA")), TOL_F32)
  expect_lte(rt_max_abs_diff(reducedDim(sce, "UMAP"), reducedDim(back, "UMAP")), TOL_F32)
})

test_that("complex SCE round-trips colData (factor, NA, character) and rowData", {
  skip_if_no_anndata()
  skip_if_not_installed("BiocNeighbors")
  sce  <- build_complex_sce()
  back <- rt_sce(sce)

  cd0 <- colData(sce); cd1 <- colData(back)
  expect_equal(as.character(cd1$cell_type), as.character(cd0$cell_type))
  expect_true(is.na(cd1$qc[1]))
  expect_equal(as.character(cd1$batch), as.character(cd0$batch))

  rd0 <- rowData(sce); rd1 <- rowData(back)
  expect_equal(as.character(rd1$symbol), as.character(rd0$symbol))
  expect_equal(as.logical(rd1$is_hvg), rd0$is_hvg)

  # Cell / feature names survive.
  expect_equal(colnames(back), colnames(sce))
  expect_equal(rownames(back), rownames(sce))
})

# ----------------------------------------------------------------------
# SCE -> AnnData: multimodal / graph data carried one direction
# ----------------------------------------------------------------------

test_that("SCE colPair becomes obsp on the AnnData side", {
  skip_if_no_anndata()
  skip_if_not_installed("BiocNeighbors")
  sce <- build_complex_sce()
  ad  <- rt_quiet(convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE))
  # extract_pairs turns the single-mcol SelfHits into a weighted sparse matrix.
  expect_true("knn" %in% rt_keys(ad$obsp))
  expect_equal(dim(as.matrix(ad$obsp[["knn"]])), c(ncol(sce), ncol(sce)))
})

test_that("SCE altExp is preserved into AnnData uns$altExperiments", {
  skip_if_no_anndata()
  skip_if_not_installed("BiocNeighbors")
  sce <- build_complex_sce()
  ad  <- rt_quiet(convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE))

  expect_true("altExperiments" %in% rt_keys(ad$uns))
  expect_true("ADT" %in% names(ad$uns[["altExperiments"]]))
  alt <- ad$uns[["altExperiments"]][["ADT"]]
  # The ADT altExp is itself an AnnData (4 ADT features x n_cells).
  expect_equal(as.integer(alt$n_vars), 4L)
  expect_equal(as.integer(alt$n_obs), ncol(sce))
})

test_that("useAltExp = FALSE drops altExps instead of storing them", {
  skip_if_no_anndata()
  skip_if_not_installed("BiocNeighbors")
  sce <- build_complex_sce()
  ad  <- rt_quiet(convert_to_anndata(sce, assayName = "counts", useAltExp = FALSE))
  expect_false("altExperiments" %in% rt_keys(ad$uns))
})

# ----------------------------------------------------------------------
# Known gaps (documented, not silently missing) -- see TODO(#18)
# ----------------------------------------------------------------------

test_that("[known gap] altExps are not restored as SCE altExps on the return leg", {
  # TODO(#18): SCE altExp -> AnnData uns$altExperiments works, but
  # convert_anndata_to_sce() does not read uns$altExperiments back, so altExps
  # do not survive SCE -> AnnData -> SCE. A reverse altExp restore step would
  # close this.
  skip("Known gap: convert_anndata_to_sce() does not restore uns$altExperiments.")
  sce  <- build_complex_sce()
  back <- rt_sce(sce)
  expect_true("ADT" %in% altExpNames(back))
})

test_that("[known gap] colPairs are not restored on the return leg", {
  # TODO(#18): SCE colPair -> AnnData obsp works, but convert_anndata_to_sce()
  # has no obsp -> colPairs step, so colPairs do not survive
  # SCE -> AnnData -> SCE (mirror of the AnnData-side obsp gap).
  skip("Known gap: convert_anndata_to_sce() has no obsp -> colPairs step.")
  sce  <- build_complex_sce()
  back <- rt_sce(sce)
  expect_gt(length(colPairNames(back)), 0)
})
