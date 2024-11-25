#' Convert Seurat Graph to colPair-Compatible Format
#'
#' This function converts a graph (adjacency matrix or edge list) into a format
#' compatible with the `colPair` slot in a `SingleCellExperiment` object. The graph
#' must match the cell indices or names in the `SingleCellExperiment` object.
#'
#' @param graph A graph represented as either a square adjacency matrix, distance matrix, or an edge list (data frame or two-column matrix).
#' @param colnames_sce A character vector of cell names from the `colData` of the `SingleCellExperiment` object.
#'                      Used to map the graph indices to the correct cell names.
#' @return A `SelfHits` object compatible with the `colPair` slot in `SingleCellExperiment`.
#' @importFrom S4Vectors SelfHits
#' @import Matrix
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
convert_graph_to_colPair <- function(graph, colnames_sce) {
  if (is.matrix(graph) || inherits(graph, "Matrix")) {
    # Ensure the matrix is square
    if (nrow(graph) != ncol(graph)) {
      stop("The graph matrix must be square.")
    }

    # Convert adjacency matrix to an edge list
    edge_list <- which(graph != 0, arr.ind = TRUE)
    from <- rownames(graph)[edge_list[, 1]]
    to <- colnames(graph)[edge_list[, 2]]
  } else if (is.data.frame(graph) || is.matrix(graph)) {
    # Assume graph is an edge list
    from <- graph[, 1]
    to <- graph[, 2]
  } else {
    stop("Unsupported graph format. Provide a square matrix or an edge list.")
  }

  # Convert cell names to indices
  from <- match(from, colnames_sce)
  to <- match(to, colnames_sce)

  # Validate indices
  if (any(is.na(from)) || any(is.na(to))) {
    stop("Graph indices must match column names in the SCE object.")
  }

  # Create and return a SelfHits object
  SelfHits(from = from, to = to, nnode = length(colnames_sce))
}
