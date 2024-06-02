#' Command Line Interface for convert2anndata
#'
#' This function serves as a command line interface for the convert2anndata package.
#' It parses command line arguments and calls the appropriate functions to convert
#' a SingleCellExperiment or Seurat object to an AnnData object.
#'
#' @import optparse
#' @importFrom anndata write_h5ad
#' @importFrom optparse make_option OptionParser parse_args
#' @export
cli_convert <- function() {
  anndata::write_h5ad

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
      call. = FALSE
    )
  }

  # Set output filename
  if (is.null(opt$output)) {
    opt$output <- sub("\\.[rR][dD][sS]$", ".h5ad", opt$input, ignore.case = TRUE)
  }

  # Load data
  timestamped_cat("Loading data from:", opt$input, "\n")
  data <- readRDS(opt$input)

  # Convert to SingleCellExperiment if necessary
  sce <- convert_seurat_to_sce(data)
  timestamped_cat("Data loaded and converted successfully if needed.\n")

  # Convert SCE to AnnData
  ad <- convert_to_anndata(sce, opt$assay, useAltExp = !isTRUE(opt$`disable-recursive-altExp`))

  # Save AnnData object
  timestamped_cat("Saving the AnnData object to:", opt$output, "\n")
  write_h5ad(ad, opt$output)
  timestamped_cat("Conversion complete:", opt$output, "\n")
}
