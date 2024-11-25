library(testthat)
library(convert2anndata)
library(SingleCellExperiment)
library(Seurat)

# Test conversion with a Seurat v5 object without layers
test_that("convert_seurat_to_sce works with Seurat v5 object without layers", {
  skip_if_not_installed("Seurat")
  # Create a Seurat v5 object
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)
  # Set version to Seurat v5
  seurat_obj@version <- package_version("5.0.0")
  sce <- convert_seurat_to_sce(seurat_obj)
  expect_s4_class(sce, "SingleCellExperiment")
})

# Test conversion with a Seurat v5 object with layers
test_that("convert_seurat_to_sce works with Seurat v5 object with layers", {
  skip_if_not_installed("Seurat")

  # Create a counts matrix
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)

  # Create an Assay5 object with layers
  assay <- CreateAssay5Object(counts = counts_matrix)

  # Define layer names and split assay into layers
  layer_names <- c(rep("Layer1", 2), rep("Layer2", 2))
  assay <- split(assay, f = layer_names)

  # Add the assay to a Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)
  seurat_obj@version <- package_version("5.0.0")
  seurat_obj@assays$RNA <- assay

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")
  expect_true("RNA_counts" %in% assayNames(sce))

  # Check that the converted layers are combined correctly
  expect_equal(as.matrix(assay(sce, "RNA_counts")), counts_matrix)

  # Validate colData contains the correct layer names
  expect_true("RNA_counts_split" %in% colnames(colData(sce)))
  expect_equal(colData(sce)$RNA_counts_split, c(rep("Layer1", 2), rep("Layer2", 2)))
})

# Test that dimensional reductions are preserved
test_that("convert_seurat_to_sce preserves dimensional reductions", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object with PCA and UMAP
  counts_matrix <- matrix(rpois(500, lambda = 10), nrow = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:50)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)

  # Preprocess data
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)

  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = 5)

  # Dynamically set n_neighbors to avoid exceeding dataset size
  n_cells <- ncol(seurat_obj)
  if (n_cells > 1) {
    n_neighbors <- min(5, n_cells - 1) # Ensure n_neighbors < number of cells

    # Run UMAP
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:5, n_neighbors = n_neighbors)
  } else {
    warning("Not enough cells to run UMAP. Skipping UMAP step.")
  }

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")

  # Check that reduced dimensions are transferred
  expect_true("pca" %in% reducedDimNames(sce))
  if (n_cells > 1) {
    expect_true("umap" %in% reducedDimNames(sce))
  } else {
    expect_false("umap" %in% reducedDimNames(sce))
  }
})

# Test that graphs (e.g., distance matrices) are transferred
test_that("convert_seurat_to_sce transfers graphs (distance matrices)", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object with a graph
  counts_matrix <- matrix(rpois(500, lambda = 10), nrow = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:50)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)

  # Preprocess data
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 5)

  # Find neighbors to create a graph
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:5)

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")

  # Check that graphs are transferred
  expect_true(length(colPairNames(sce)) > 0)
  expect_true("RNA_nn" %in% colPairNames(sce) | "SNN" %in% colPairNames(sce))

  # Optionally check the content of the graph
  graph_name <- colPairNames(sce)[1]
  expect_true(!is.null(colPair(sce, graph_name)))
})

# Test that misc and tools data are transferred
test_that("convert_seurat_to_sce transfers misc and tools data", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object
  counts_matrix <- matrix(rpois(500, lambda = 10), nrow = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:50)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)

  # Add data to misc and tools slots
  seurat_obj@misc$custom_data <- list(a = 1, b = 2)
  seurat_obj@tools$analysis_result <- data.frame(result = rnorm(50))

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")

  # Check that misc data is transferred
  expect_true("misc" %in% names(metadata(sce)))
  expect_equal(metadata(sce)$misc$custom_data, seurat_obj@misc$custom_data)

  # Check that tools data is transferred
  expect_true("tools" %in% names(metadata(sce)))
  expect_equal(metadata(sce)$tools$analysis_result, seurat_obj@tools$analysis_result)
})

# Test that commands history and version information are transferred
test_that("convert_seurat_to_sce transfers commands history information", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object
  counts_matrix <- matrix(rpois(500, lambda = 10), nrow = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:50)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)

  # Run some commands to populate the commands slot
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")

  # Check that commands history is transferred
  expect_true("commands" %in% names(metadata(sce)))
  expect_equal(metadata(sce)$commands, seurat_obj@commands)
})

# Test that additional unstructured data is transferred
test_that("convert_seurat_to_sce transfers additional unstructured data", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object
  counts_matrix <- matrix(rpois(500, lambda = 10), nrow = 10)
  rownames(counts_matrix) <- paste0("Gene", 1:10)
  colnames(counts_matrix) <- paste0("Cell", 1:50)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)

  # Add custom data to the @misc slot (simulate unstructured data)
  seurat_obj@misc$custom_data <- list(custom_field = "custom_value")

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")

  # Check that the custom data in @misc is transferred to metadata
  expect_true("misc" %in% names(metadata(sce)))
  expect_equal(metadata(sce)$misc$custom_data, seurat_obj@misc$custom_data)
})

# Test conversion when input is already a SingleCellExperiment
test_that("convert_seurat_to_sce works with SingleCellExperiment objects", {
  sce_input <- SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
  sce_output <- convert_seurat_to_sce(sce_input)
  expect_s4_class(sce_output, "SingleCellExperiment")
  # Check that the output is identical to the input
  expect_equal(sce_input, sce_output)
})

# Test handling of unsupported object types
test_that("convert_seurat_to_sce handles unsupported object types gracefully", {
  unsupported_obj <- data.frame(a = 1:5, b = 6:10)
  expect_error(convert_seurat_to_sce(unsupported_obj), "Failed to convert object to SingleCellExperiment")
})

# Test handling of NULL input
test_that("convert_seurat_to_sce handles NULL input gracefully", {
  expect_error(convert_seurat_to_sce(NULL), "Failed to convert object to SingleCellExperiment")
})

# Test handling of Seurat object with mismatching assays
test_that("convert_seurat_to_sce handles Seurat object with mismatching assays", {
  skip_if_not_installed("Seurat")
  # Create a Seurat object with additional assays having different features
  counts_main <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_main) <- paste0("Gene", 1:5)
  colnames(counts_main) <- paste0("Cell", 1:4)
  seurat_obj <- CreateSeuratObject(counts = counts_main)
  counts_assay2 <- matrix(1:16, nrow = 4, ncol = 4)
  rownames(counts_assay2) <- paste0("Gene", 6:9) # Different genes
  colnames(counts_assay2) <- colnames(counts_main)
  seurat_obj[["Assay2"]] <- CreateAssayObject(counts = counts_assay2)
  sce <- convert_seurat_to_sce(seurat_obj)
  expect_s4_class(sce, "SingleCellExperiment")
  # Check that altExps contain the mismatching assay
  expect_true("Assay2" %in% altExpNames(sce))
})

# Test handling of Seurat object with layers but minimal data
test_that("convert_seurat_to_sce handles Seurat object with layers but minimal data", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)

  # Simulate Seurat v5 object
  seurat_obj@version <- package_version("5.0.0")

  # Add a minimal data layer
  minimal_layer <- matrix(0, nrow = 5, ncol = 4)
  rownames(minimal_layer) <- rownames(seurat_obj)
  colnames(minimal_layer) <- colnames(seurat_obj)
  seurat_obj@assays$RNA@layers$counts <- as(minimal_layer, "dgCMatrix") # Layer name should be "counts"

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")
  # Since the layer has minimal data, it may be skipped
})

# Test handling of Seurat object with multiple assays
test_that("convert_seurat_to_sce handles Seurat object with multiple assays", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object with multiple assays
  counts_main <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_main) <- paste0("Gene", 1:5)
  colnames(counts_main) <- paste0("Cell", 1:4)
  seurat_obj <- CreateSeuratObject(counts = counts_main)

  # Add an additional assay with a different feature set
  counts_assay2 <- matrix(1:16, nrow = 4, ncol = 4)
  rownames(counts_assay2) <- paste0("Gene", 6:9) # Different genes
  colnames(counts_assay2) <- colnames(counts_main)
  seurat_obj[["Assay2"]] <- CreateAssayObject(counts = counts_assay2)

  # Add another assay with yet another feature set
  counts_assay3 <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(counts_assay3) <- paste0("Gene", 10:12) # Different genes
  colnames(counts_assay3) <- colnames(counts_main)
  seurat_obj[["Assay3"]] <- CreateAssayObject(counts = counts_assay3)

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")

  # Check that altExps contain the mismatching assays
  expect_true("Assay2" %in% altExpNames(sce))
  expect_true("Assay3" %in% altExpNames(sce))
})

# Test handling of Seurat object with per-cell metadata
test_that("convert_seurat_to_sce handles Seurat object with per-cell metadata", {
  skip_if_not_installed("Seurat")

  # Create a Seurat object with per-cell metadata
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)
  seurat_obj$cell_type <- c("A", "B", "A", "B")

  # Convert to SingleCellExperiment
  sce <- convert_seurat_to_sce(seurat_obj)

  # Validate output
  expect_s4_class(sce, "SingleCellExperiment")

  # Check that colData contains the metadata
  expect_true("cell_type" %in% colnames(colData(sce)))

  # Compare metadata values
  expect_equal(
    unname(colData(sce)$cell_type),
    unname(seurat_obj$cell_type) # Unset names for comparison
  )
})


# Test handling of invalid input types
test_that("convert_seurat_to_sce handles invalid input types", {
  expect_error(convert_seurat_to_sce(123), "Failed to convert object to SingleCellExperiment")
  expect_error(convert_seurat_to_sce("invalid_input"), "Failed to convert object to SingleCellExperiment")
})



# Test handling of Seurat object with missing metadata
test_that("convert_seurat_to_sce handles Seurat object with missing metadata", {
  skip_if_not_installed("Seurat")
  # Create a Seurat object without metadata
  counts_matrix <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:4)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)
  # Remove meta.data
  seurat_obj@meta.data <- data.frame(row.names = colnames(seurat_obj))
  sce <- convert_seurat_to_sce(seurat_obj)
  expect_s4_class(sce, "SingleCellExperiment")
  # Check that colData contains only default columns
  expect_true(ncol(colData(sce)) == 0 || all(colnames(colData(sce)) %in% c("orig.ident")))
})
