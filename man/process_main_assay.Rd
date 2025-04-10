% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_main_assay.R
\name{process_main_assay}
\alias{process_main_assay}
\title{Process Main Assay}
\usage{
process_main_assay(sce, assayName)
}
\arguments{
\item{sce}{A SingleCellExperiment object.}

\item{assayName}{The name of the assay to use as the main data matrix in AnnData.}
}
\value{
The main data matrix (anndata.X) for the AnnData object.
}
\description{
This function processes the main assay in the SingleCellExperiment object to be used as the primary data matrix in the AnnData object.
}
