% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_graph_to_colPair.R
\name{convert_graph_to_colPair}
\alias{convert_graph_to_colPair}
\title{Convert Seurat Graph to colPair-Compatible Format}
\usage{
convert_graph_to_colPair(graph, colnames_sce)
}
\arguments{
\item{graph}{A graph represented as either a square adjacency matrix, distance matrix, or an edge list (data frame or two-column matrix).}

\item{colnames_sce}{A character vector of cell names from the `colData` of the `SingleCellExperiment` object.
Used to map the graph indices to the correct cell names.}
}
\value{
A `SelfHits` object compatible with the `colPair` slot in `SingleCellExperiment`.
}
\description{
This function converts a graph (adjacency matrix or edge list) into a format
compatible with the `colPair` slot in a `SingleCellExperiment` object. The graph
must match the cell indices or names in the `SingleCellExperiment` object.
}
