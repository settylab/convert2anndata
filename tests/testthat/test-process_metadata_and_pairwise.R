library(testthat)
library(SingleCellExperiment)
library(convert2anndata)

test_that("process_metadata_and_pairwise handles metadata correctly", {
  # Create a basic SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:6, ncol = 2)))
  X <- matrix(1:6, ncol = 2)
  n_obs <- nrow(X)
  n_var <- ncol(X)
  
  # Add metadata
  metadata(sce)$obs_pairwise <- matrix(1, nrow = n_obs, ncol = n_obs)
  metadata(sce)$var_pairwise <- matrix(1, nrow = n_var, ncol = n_var)
  metadata(sce)$obs_matrix <- matrix(1, nrow = n_obs, ncol = 5)
  metadata(sce)$var_matrix <- matrix(1, nrow = n_var, ncol = 5)
  
  # Process metadata
  uns_data <- process_metadata_and_pairwise(sce, list(), X)
  
  # Check outputs
  expect_true(is.list(uns_data$uns))
  expect_true(is.list(uns_data$obsp))
  expect_true(is.list(uns_data$varp))
  expect_true(is.list(uns_data$varm))
  expect_true(is.list(uns_data$obsm))
  
  # Validate transferred metadata
  expect_true("obs_pairwise" %in% names(uns_data$obsp))
  expect_true("var_pairwise" %in% names(uns_data$varp))
  expect_true("obs_matrix" %in% names(uns_data$obsm))
  expect_true("var_matrix" %in% names(uns_data$varm))
  
  # Validate remaining metadata
  expect_false("obs_pairwise" %in% names(uns_data$uns))
  expect_false("var_pairwise" %in% names(uns_data$uns))
  expect_false("obs_matrix" %in% names(uns_data$uns))
  expect_false("var_matrix" %in% names(uns_data$uns))
})

test_that("process_metadata_and_pairwise handles alternative experiments", {
  # Create a basic SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  X <- matrix(1:4, ncol = 2)
  
  # Add alternative experiments
  alt_exps <- list(alt1 = list(data = "altExp1"), alt2 = list(data = "altExp2"))
  
  # Process metadata with alternative experiments
  uns_data <- process_metadata_and_pairwise(sce, alt_exps, X)
  
  # Check altExperiments in uns
  expect_true("altExperiments" %in% names(uns_data$uns))
  expect_equal(uns_data$uns$altExperiments, alt_exps)
})

test_that("process_metadata_and_pairwise handles empty or missing metadata", {
  # Create a SingleCellExperiment without metadata
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  X <- matrix(1:4, ncol = 2)
  
  # Process metadata
  uns_data <- process_metadata_and_pairwise(sce, list(), X)
  
  # Check outputs
  expect_true(is.list(uns_data$uns))
  expect_equal(names(uns_data$uns), c("version"))
  expect_equal(length(uns_data$obsp), 0)
  expect_equal(length(uns_data$varp), 0)
  expect_equal(length(uns_data$varm), 0)
  expect_equal(length(uns_data$obsm), 0)
})

test_that("process_metadata_and_pairwise handles mismatched dimensions", {
  # Create a SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  X <- matrix(1:4, ncol = 2)
  n_obs <- nrow(X)
  n_var <- ncol(X)
  
  # Add mismatched metadata
  metadata(sce)$mismatched <- matrix(1, nrow = n_obs + 1, ncol = n_obs + 1)
  
  # Process metadata
  uns_data <- process_metadata_and_pairwise(sce, list(), X)
  
  # Check mismatched metadata is not transferred
  expect_true("mismatched" %in% names(uns_data$uns))
  expect_false("mismatched" %in% names(uns_data$obsp))
  expect_false("mismatched" %in% names(uns_data$varp))
  expect_false("mismatched" %in% names(uns_data$obsm))
  expect_false("mismatched" %in% names(uns_data$varm))
})