#!/usr/bin/env Rscript

# Load helper functions
source(system.file("R/utils.R", package = "yourpackage"))

# List of required packages
required_packages <- c("SingleCellExperiment", "anndata", "optparse")
check_and_load_packages(required_packages)

# Description and help
description <- paste(
  "This script converts a potentially old Seurat object or a SingleCellExperiment",
  "stored in an RDS file into an AnnData object stored as an H5AD file.",
  "The user can specify input and output file paths, with an option to change",
  "the output filename from .rds to .h5ad if no output is specified."
)

# Set up command-line options
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Path to the input RDS file containing the SingleCellExperiment object. This option is required.",
    metavar = "file"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = paste(
      "Path to the output H5AD file. If not specified,",
      "the output path is derived by replacing the .rds",
      "extension of the input path with .h5ad."
    ),
    metavar = "file"
  ),
  make_option(c("-a", "--assay"),
    type = "character", default = "counts",
    help = "The assay to use as the main matrix (anndata.X). Defaults to 'counts'.",
    metavar = "assay_name"
  ),
  make_option(c("-d", "--disable-recursive-altExp"),
    action = "store_true", default = FALSE,
    help = "Disable recursive recovery of altExperiments and discard them instead.",
    metavar = "boolean"
  )
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list, description = description)
opt <- parse_args(opt_parser)

# Check if input file is provided
if (is.null(opt$input)) {
  stop("No input file provided. Use --input to specify the RDS file.", 
       call. = FALSE)
}

# Set output filename
if (is.null(opt$output)) {
  opt$output <- sub("\\.[rR][dD][sS]$", ".h5ad", opt$input, ignore.case = TRUE)
}

# Function to convert SingleCellExperiment to AnnData
convert_to_anndata <- function(sce, assay_name) {
  # Function implementation here...
}

timestamped_cat("Loading data from:", opt$input, "\n")
data <- readRDS(opt$input)

data <- convert_seurat_to_sce(data)
timestamped_cat("Data loaded and converted successfully if needed.\n")

# Convert SCE to AnnData
ad <- convert_to_anndata(sce, opt$assay, useAltExp=!isTRUE(opt$`disable-recursive-altExp`))

timestamped_cat("Saving the AnnData object to:", opt$output, "\n")
write_h5ad(ad, opt$output)
timestamped_cat("Conversion complete: ", opt$output, "\n")
