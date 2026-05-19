library(testthat)
library(convert2anndata)
library(Matrix)
library(anndata)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)

# ======================================================================
# Complex AnnData round-trip fidelity (issue #18).
#
# Two full cycles are exercised on a COMPLEX synthetic AnnData (two layers,
# PCA + UMAP obsm, varm, two obsp graphs, varp, raw, and obs/var carrying
# factor / NA / logical / character / numeric columns):
#
#   Path A  AnnData -> SCE   -> AnnData            (SCE-mediated, faithful)
#   Path B  AnnData -> Seurat -> SCE -> AnnData    (Seurat-mediated, lossy)
#
# Fixtures and helpers live in helper-roundtrip.R. Path A writes through an
# on-disk .h5ad so serialisation is part of the test. Components that the
# converters do not (yet) carry in a given direction are documented as
# skip()'d tests with a TODO(#18) so the gap is recorded, not silent.
# ======================================================================

TOL     <- 1e-6   # exact-ish for values that stay float64 across the trip
TOL_F32 <- 1e-4   # looser bound where a value passes through float32

# ----------------------------------------------------------------------
# Path A: AnnData -> SCE -> AnnData  (the faithful, SCE-mediated path)
# ----------------------------------------------------------------------

test_that("Path A: complex AnnData round-trips X/counts through SCE and h5ad", {
  skip_if_no_anndata()
  ad <- build_complex_anndata()
  sce <- rt_quiet(convert_anndata_to_sce(ad, counts_layer = "counts"))
  ad_rt <- rt_write_read(rt_quiet(
    convert_to_anndata(sce, assayName = "counts", useAltExp = TRUE)
  ))

  # The original counts layer should land in the round-tripped X.
  expect_equal(rt_max_abs_diff(ad$layers[["counts"]], ad_rt$X), 0)
  # obs/var names survive the full trip.
  expect_identical(as.character(ad$obs_names), as.character(ad_rt$obs_names))
  expect_identical(as.character(ad$var_names), as.character(ad_rt$var_names))
})

test_that("Path A: extra layers survive by name (counts + logcounts)", {
  skip_if_no_anndata()
  ad <- build_complex_anndata()
  sce <- rt_quiet(convert_anndata_to_sce(ad, counts_layer = "counts"))
  ad_rt <- rt_write_read(rt_quiet(convert_to_anndata(sce, assayName = "counts")))

  expect_true("logcounts" %in% rt_keys(ad_rt$layers))
  expect_lte(rt_max_abs_diff(ad$layers[["logcounts"]],
                             ad_rt$layers[["logcounts"]]), TOL)
})

test_that("Path A: obsm reductions (PCA + UMAP) survive within tolerance", {
  skip_if_no_anndata()
  ad <- build_complex_anndata()
  sce <- rt_quiet(convert_anndata_to_sce(ad, counts_layer = "counts"))
  ad_rt <- rt_write_read(rt_quiet(convert_to_anndata(sce, assayName = "counts")))

  rt_map <- setNames(rt_keys(ad_rt$obsm), rt_norm_key(rt_keys(ad_rt$obsm)))
  for (k in c("X_pca", "X_umap")) {
    nk <- rt_norm_key(k)
    expect_true(nk %in% names(rt_map),
                info = sprintf("obsm '%s' missing after round trip", k))
    expect_lte(rt_max_abs_diff(ad$obsm[[k]], ad_rt$obsm[[rt_map[[nk]]]]), TOL_F32)
  }
})

test_that("Path A: obs metadata (factor, NA, logical, character) survives", {
  skip_if_no_anndata()
  ad <- build_complex_anndata()
  sce <- rt_quiet(convert_anndata_to_sce(ad, counts_layer = "counts"))
  ad_rt <- rt_write_read(rt_quiet(convert_to_anndata(sce, assayName = "counts")))

  obs0 <- as.data.frame(ad$obs)
  obs1 <- as.data.frame(ad_rt$obs)

  # Categorical preserved (levels + per-cell labels), NA kept as NA.
  expect_equal(as.character(obs1$cell_type), as.character(obs0$cell_type))
  expect_true(is.na(obs1$pct_mt[1]))
  expect_equal(obs1$n_counts, obs0$n_counts)
  expect_equal(as.logical(obs1$is_ctrl), obs0$is_ctrl)
  expect_equal(as.character(obs1$donor), obs0$donor)
})

test_that("Path A: var metadata (factor, logical, numeric) survives", {
  skip_if_no_anndata()
  ad <- build_complex_anndata()
  sce <- rt_quiet(convert_anndata_to_sce(ad, counts_layer = "counts"))
  ad_rt <- rt_write_read(rt_quiet(convert_to_anndata(sce, assayName = "counts")))

  var0 <- as.data.frame(ad$var)
  var1 <- as.data.frame(ad_rt$var)
  expect_equal(as.character(var1$gene_class), as.character(var0$gene_class))
  expect_equal(as.logical(var1$highly_variable), var0$highly_variable)
  expect_equal(var1$mean_expr, var0$mean_expr, tolerance = TOL_F32)
})

test_that("Path A: both sparse and dense X inputs round-trip counts faithfully", {
  skip_if_no_anndata()
  for (sparse in c(TRUE, FALSE)) {
    ad <- build_complex_anndata(sparse_X = sparse)
    sce <- rt_quiet(convert_anndata_to_sce(ad, counts_layer = "counts"))
    ad_rt <- rt_write_read(rt_quiet(convert_to_anndata(sce, assayName = "counts")))
    expect_equal(rt_max_abs_diff(ad$layers[["counts"]], ad_rt$X), 0,
                 info = sprintf("sparse_X = %s", sparse))
  }
})

# ----------------------------------------------------------------------
# AnnData -> Seurat: complex components attach directly
# ----------------------------------------------------------------------

test_that("complex AnnData attaches reductions, obsp graphs, and var feature-meta to Seurat", {
  skip_if_no_anndata()
  ad <- build_complex_anndata()
  s <- rt_quiet(convert_anndata_to_seurat(ad, counts_layer = "counts"))

  expect_setequal(names(s@reductions), c("pca", "umap"))
  expect_true(all(c("RNA_connectivities", "RNA_distances") %in% names(s@graphs)))
  fm <- s[["RNA"]][[]]
  expect_true(all(c("gene_class", "highly_variable", "mean_expr") %in% colnames(fm)))
  # PCA embeddings preserved on the direct AnnData -> Seurat hop.
  expect_lte(rt_max_abs_diff(s@reductions$pca@cell.embeddings, ad$obsm[["X_pca"]]),
             TOL_F32)
})

test_that("Path B: obsp connectivity survives the Seurat-mediated round trip", {
  skip_if_no_anndata()
  # AnnData obsp -> Seurat Graph -> SCE colPair -> AnnData obsp. Edge *weights*
  # are dropped by convert_graph_to_colPair (connectivity is binarised), so we
  # assert the non-zero PATTERN is preserved rather than the weights.
  ad <- build_complex_anndata()
  s   <- rt_quiet(convert_anndata_to_seurat(ad, counts_layer = "counts"))
  sce <- rt_quiet(convert_seurat_to_sce(s))
  ad_b <- rt_quiet(convert_to_anndata(sce, assayName = "counts"))

  # The Seurat hop prefixes the graph with the assay name
  # (connectivities -> RNA_connectivities), so look it up by the prefixed key.
  expect_true("RNA_connectivities" %in% rt_keys(ad_b$obsp))
  orig <- as.matrix(ad$obsp[["connectivities"]])
  back <- as.matrix(ad_b$obsp[["RNA_connectivities"]])
  expect_equal(dim(back), dim(orig))
  expect_equal(rt_pattern(back), rt_pattern(orig))
})

# ----------------------------------------------------------------------
# Known gaps (documented, not silently missing) -- see TODO(#18)
# ----------------------------------------------------------------------

test_that("[known gap] Path B drops extra AnnData layers (Seurat assay model)", {
  # TODO(#18): a Seurat assay has only counts/data/scale.data slots, so an
  # extra layer such as 'logcounts' has nowhere to live on the Seurat hop.
  # Use the SCE-mediated path (Path A) when verbatim layers must survive.
  skip("Known structural gap: Seurat-mediated path cannot carry extra layers.")
  ad <- build_complex_anndata()
  s   <- convert_anndata_to_seurat(ad, counts_layer = "counts")
  sce <- convert_seurat_to_sce(s)
  ad_b <- convert_to_anndata(sce, assayName = "counts")
  expect_true("logcounts" %in% rt_keys(ad_b$layers))
})

test_that("[known gap] Path B drops var feature metadata", {
  # TODO(#18): convert_anndata_to_seurat() attaches var to the Seurat feature
  # meta, but convert_seurat_to_sce() does not lift feature meta back into
  # rowData, so 'var' is lost across AnnData -> Seurat -> SCE -> AnnData.
  # Use the SCE-mediated path (Path A) when var metadata must survive.
  skip("Known structural gap: Seurat-mediated path does not preserve var.")
  ad <- build_complex_anndata()
  s   <- convert_anndata_to_seurat(ad, counts_layer = "counts")
  sce <- convert_seurat_to_sce(s)
  ad_b <- convert_to_anndata(sce, assayName = "counts")
  expect_true("gene_class" %in% rt_keys(ad_b$var))
})

test_that("[known gap] convert_anndata_to_sce does not wire obsp into colPairs", {
  # TODO(#18): convert_anndata_to_seurat() attaches obsp as Graphs and
  # extract_anndata_obsp() already returns the matrices, but
  # convert_anndata_to_sce() has no obsp -> colPairs step, so the SCE-mediated
  # AnnData round trip loses obsp. A weight-preserving SelfHits construction
  # would close this (see the issue report).
  skip("Known gap: convert_anndata_to_sce() has no obsp -> colPairs step.")
  ad <- build_complex_anndata()
  sce <- convert_anndata_to_sce(ad, counts_layer = "counts")
  expect_gt(length(SingleCellExperiment::colPairNames(sce)), 0)
})

test_that("[known gap] convert_anndata_to_sce does not carry varm / varp / raw", {
  # TODO(#18): convert_anndata_to_sce() reads only X/layers/obs/var/obsm, so
  # varm (gene loadings), varp (feature graphs) and adata.raw are dropped on
  # the AnnData -> SCE hop and therefore do not survive the SCE-mediated round
  # trip. These are unimplemented reverse directions, not regressions.
  skip("Known gap: convert_anndata_to_sce() ignores varm/varp/raw.")
  ad <- build_complex_anndata()
  sce <- convert_anndata_to_sce(ad, counts_layer = "counts")
  ad_rt <- convert_to_anndata(sce, assayName = "counts")
  expect_true("PCs" %in% rt_keys(ad_rt$varm))
  expect_true("corr" %in% rt_keys(ad_rt$varp))
})
