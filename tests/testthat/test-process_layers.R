library(testthat)
library(Seurat)

test_that("process_layers combines split layers with collisions correctly", {
  skip_if_not_installed("Seurat")

  # Create a counts matrix
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)

  # Create an Assay5 object with split layers
  assay <- CreateAssay5Object(counts = counts_matrix)
  layer1 <- counts_matrix[, 1:2]
  layer2 <- counts_matrix[, 3:4]
  assay <- SetAssayData(assay, layer = "counts.A", new.data = layer1)
  assay <- SetAssayData(assay, layer = "counts.B", new.data = layer2)

  # Process layers
  result <- process_layers(assay)

  # Validate split names
  expect_equal(result$split_names[["counts"]], c("A", "A", "B", "B"))

  # Validate combined counts
  combined_counts <- GetAssayData(result$assay_object, layer = "counts.1")
  expect_equal(as.matrix(combined_counts), counts_matrix)

  # Check that the original split layers are removed
  remaining_layers <- Layers(result$assay_object)
  expect_false("counts.A" %in% remaining_layers)
  expect_false("counts.B" %in% remaining_layers)
  expect_true("counts.1" %in% remaining_layers)
  expect_true("counts" %in% remaining_layers)
})

test_that("process_layers combines split layers without collisions correctly", {
  skip_if_not_installed("Seurat")

  # Create a counts matrix
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)

  # Create an Assay5 object with split layers
  assay <- CreateAssay5Object(counts = counts_matrix)
  splitting <- c("A", "B", "A", "B")
  assay <- split(assay, f = c("A", "B", "A", "B"))

  # Process layers
  result <- process_layers(assay)

  # Validate split names
  expect_equal(result$split_names[["counts"]], splitting)

  # Validate combined counts
  combined_counts <- GetAssayData(result$assay_object, layer = "counts")
  expect_equal(as.matrix(combined_counts), counts_matrix)

  # Check that the original split layers are removed
  remaining_layers <- Layers(result$assay_object)
  expect_false("counts.A" %in% remaining_layers)
  expect_false("counts.B" %in% remaining_layers)
  expect_false("counts.1" %in% remaining_layers)
  expect_true("counts" %in% remaining_layers)
})

test_that("process_layers handles non-split layers correctly", {
  skip_if_not_installed("Seurat")

  # Create a counts matrix
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)

  # Create an Assay5 object with a single layer
  assay <- CreateAssay5Object(counts = counts_matrix)

  # Add a single non-split layer
  assay <- SetAssayData(assay, layer = "data", new.data = counts_matrix)

  # Process layers
  result <- process_layers(assay)

  # Validate that split_names is an empty list
  expect_equal(result$split_names, list())

  # Validate that the non-split layer remains intact
  remaining_layers <- Layers(result$assay_object)
  expect_true("data" %in% remaining_layers)
  expect_equal(length(remaining_layers), 2)
})

test_that("process_layers handles mixed split and non-split layers", {
  skip_if_not_installed("Seurat")

  # Create a counts matrix
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)

  # Create an Assay5 object with split and non-split layers
  assay <- CreateAssay5Object(counts = counts_matrix)
  layer1 <- counts_matrix[, 1:2]
  layer2 <- counts_matrix[, 3:4]
  assay <- SetAssayData(assay, layer = "test.A", new.data = layer1)
  assay <- SetAssayData(assay, layer = "test.B", new.data = layer2)

  # Add a non-split layer
  assay <- SetAssayData(assay, layer = "scale.data", new.data = counts_matrix)

  # Process layers
  result <- process_layers(assay)

  # Validate split names for combined split layers
  expect_equal(result$split_names, list(test = c("A", "A", "B", "B")))

  # Validate combined counts
  combined_counts <- GetAssayData(result$assay_object, layer = "counts")
  expect_equal(as.matrix(combined_counts), counts_matrix)

  # Validate that the non-split layer remains intact
  remaining_layers <- Layers(result$assay_object)
  expect_true("scale.data" %in% remaining_layers)
  expect_true("counts" %in% remaining_layers)
  expect_true("test" %in% remaining_layers)
  expect_false("test.A" %in% remaining_layers)
  expect_false("test.B" %in% remaining_layers)
})
