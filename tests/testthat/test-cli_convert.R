library(testthat)
library(optparse)
library(reticulate)
library(withr)
library(convert2anndata)

# Ensure anndata is loaded
anndata <- reticulate::import("anndata", convert = FALSE)

test_that("cli_convert works with SingleCellExperiment input", {
  with_options(list(args = c("-i", "tests/testdata/mock_sce.rds", "-o", "tests/testdata/output_sce.h5ad")), {
    expect_message(cli_convert(), "Conversion complete")
    expect_true(file.exists("tests/testdata/output_sce.h5ad"))
  })
})

if (requireNamespace("Seurat", quietly = TRUE)) {
  test_that("cli_convert works with Seurat input", {
    with_options(list(args = c("-i", "tests/testdata/mock_seurat.rds", "-o", "tests/testdata/output_seurat.h5ad")), {
      expect_message(cli_convert(), "Conversion complete")
      expect_true(file.exists("tests/testdata/output_seurat.h5ad"))
    })
  })
}

test_that("cli_convert uses default assay", {
  with_options(list(args = c("-i", "tests/testdata/mock_sce.rds", "-o", "tests/testdata/output_default_assay.h5ad", "-a", "counts")), {
    expect_message(cli_convert(), "Conversion complete")
    expect_true(file.exists("tests/testdata/output_default_assay.h5ad"))
  })
})

test_that("cli_convert handles altExp flag", {
  with_options(list(args = c("-i", "tests/testdata/mock_sce.rds", "-o", "tests/testdata/output_altExp.h5ad", "-d")), {
    expect_message(cli_convert(), "Conversion complete")
    expect_true(file.exists("tests/testdata/output_altExp.h5ad"))
  })
})

test_that("cli_convert stops without input", {
  with_options(list(args = c("-o", "tests/testdata/output_no_input.h5ad")), {
    expect_error(cli_convert(), "No input file provided")
  })
})
