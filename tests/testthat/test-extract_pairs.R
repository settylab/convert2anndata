library(testthat)
library(SingleCellExperiment)
library(Matrix)
library(methods)
library(BiocNeighbors)
library(convert2anndata)

# Helper function to create a mock SingleCellExperiment object with various pairwise data scenarios
create_mock_sce_with_pairs <- function() {
  # Create a main assay matrix
  counts <- matrix(runif(100), nrow = 10, ncol = 10)

  # Create colData and rowData
  colData <- DataFrame(cell_type = rep(c("A", "B"), 5))
  rowData <- DataFrame(gene_type = rep(c("gene1", "gene2"), 5))

  # Create the main SCE object
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = colData,
    rowData = rowData
  )

  # Find KNN and add pairwise data
  knn <- findKNN(as.matrix(counts), k = 3)
  knn_index <- SelfHits(
    from = rep(1:10, each = 3), to = rep(1:3, times = 10),
    x = as.vector(knn$index), nnode = 10
  )
  knn_distance <- SelfHits(
    from = rep(1:10, each = 3), to = rep(1:3, times = 10),
    x = as.vector(knn$distance), nnode = 10
  )

  # Create a SelfHits object with multiple metadata columns
  knn_multiple_meta <- SelfHits(
    from = rep(1:10, each = 3), to = rep(1:3, times = 10),
    x = as.vector(knn$index), y = as.vector(knn$distance),
    nnode = 10
  )

  # Create a SelfHits object with no metadata columns
  knn_no_meta <- SelfHits(
    from = rep(1:10, each = 3), to = rep(1:3, times = 10),
    nnode = 10
  )

  colPairs(sce) <- list(
    index = knn_index,
    distance = knn_distance,
    knn_with_multiple_meta = knn_multiple_meta,
    knn_no_meta = knn_no_meta
  )

  return(sce)
}

test_that("extract_pairs works with various pairwise data scenarios", {
  sce <- create_mock_sce_with_pairs()

  # Extract pairs using the extract_pairs function
  pairs_list <- extract_pairs(colPairs, sce)

  # Check that the pairs_list contains the expected matrices
  expect_true("distance" %in% names(pairs_list))
  expect_true("connectivity" %in% names(pairs_list))
  expect_true("knn_with_multiple_meta_x" %in% names(pairs_list))
  expect_true("knn_with_multiple_meta_y" %in% names(pairs_list))
  expect_true("knn_no_meta" %in% names(pairs_list))

  # Verify the content of the extracted pairs
  expect_equal(dim(pairs_list[["distance"]]), c(ncol(sce), ncol(sce)))
  expect_equal(dim(pairs_list[["connectivity"]]), c(ncol(sce), ncol(sce)))
  expect_equal(dim(pairs_list[["knn_with_multiple_meta_x"]]), c(ncol(sce), ncol(sce)))
  expect_equal(dim(pairs_list[["knn_with_multiple_meta_y"]]), c(ncol(sce), ncol(sce)))
  expect_equal(dim(pairs_list[["knn_no_meta"]]), c(ncol(sce), ncol(sce)))

  # Check that the distance matrix contains the correct values
  distances <- mcols(colPairs(sce)[["distance"]])[[1]]
  indices <- mcols(colPairs(sce)[["index"]])[[1]]
  expect_equal(pairs_list[["distance"]][1, indices[1]], distances[1])
  expect_equal(pairs_list[["distance"]][1, indices[2]], distances[2])
  expect_equal(pairs_list[["distance"]][1, indices[3]], distances[3])

  # Check that the connectivity matrix contains the correct values
  expect_equal(pairs_list[["connectivity"]][1, indices[1]], 1)
  expect_equal(pairs_list[["connectivity"]][1, indices[2]], 1)
  expect_equal(pairs_list[["connectivity"]][1, indices[3]], 1)
})
