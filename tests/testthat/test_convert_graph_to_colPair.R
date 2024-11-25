# tests/testthat/test_convert_graph_to_colPair.R

library(testthat)
library(S4Vectors)
library(SingleCellExperiment)

test_that("convert_graph_to_colPair works with adjacency matrix", {
  # Create an adjacency matrix
  graph <- matrix(0, nrow = 5, ncol = 5)
  rownames(graph) <- colnames(graph) <- paste0("Cell", 1:5)
  graph[1, 2] <- 1
  graph[2, 3] <- 1
  graph[4, 5] <- 1

  # Define cell names
  colnames_sce <- paste0("Cell", 1:5)

  # Convert the graph
  result <- convert_graph_to_colPair(graph, colnames_sce)

  # Check class of the result
  expect_s4_class(result, "SelfHits")

  # Check that indices match the expected values
  expect_equal(from(result), c(1, 2, 4))
  expect_equal(to(result), c(2, 3, 5))
})

test_that("convert_graph_to_colPair works with edge list", {
  # Create an edge list
  graph <- data.frame(
    from = c("Cell1", "Cell2", "Cell4"),
    to = c("Cell2", "Cell3", "Cell5")
  )

  # Define cell names
  colnames_sce <- paste0("Cell", 1:5)

  # Convert the graph
  result <- convert_graph_to_colPair(graph, colnames_sce)

  # Check class of the result
  expect_s4_class(result, "SelfHits")

  # Check that indices match the expected values
  expect_equal(from(result), c(1, 2, 4))
  expect_equal(to(result), c(2, 3, 5))
})

test_that("convert_graph_to_colPair handles mismatched cell names", {
  # Create an edge list with mismatched names
  graph <- data.frame(
    from = c("Cell1", "CellX", "Cell4"),
    to = c("Cell2", "Cell3", "CellY")
  )

  # Define cell names
  colnames_sce <- paste0("Cell", 1:5)

  # Expect an error due to mismatched names
  expect_error(
    convert_graph_to_colPair(graph, colnames_sce),
    "Graph indices must match column names in the SCE object."
  )
})

test_that("convert_graph_to_colPair handles non-square matrix", {
  # Create a non-square matrix
  graph <- matrix(0, nrow = 4, ncol = 5)
  rownames(graph) <- paste0("Cell", 1:4)
  colnames(graph) <- paste0("Cell", 1:5)

  # Define cell names
  colnames_sce <- paste0("Cell", 1:5)

  # Expect an error due to non-square matrix
  expect_error(
    convert_graph_to_colPair(graph, colnames_sce),
    "The graph matrix must be square."
  )
})

test_that("convert_graph_to_colPair works with sparse adjacency matrix", {
  library(Matrix)

  # Create a sparse adjacency matrix
  graph <- Matrix(0, nrow = 5, ncol = 5, sparse = TRUE)
  rownames(graph) <- colnames(graph) <- paste0("Cell", 1:5)
  graph[1, 2] <- 1
  graph[3, 4] <- 1

  # Define cell names
  colnames_sce <- paste0("Cell", 1:5)

  # Convert the graph
  result <- convert_graph_to_colPair(graph, colnames_sce)

  # Check class of the result
  expect_s4_class(result, "SelfHits")

  # Check that indices match the expected values
  expect_equal(from(result), c(1, 3))
  expect_equal(to(result), c(2, 4))
})
