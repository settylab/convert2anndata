#' Convert SingleCellExperiment to AnnData
#'
#' This function converts a SingleCellExperiment (SCE) object to an AnnData object.
#' It processes assays, dimensional reductions, and metadata, including alternative experiments (altExps).
#' The main assay is used as the primary data matrix in the AnnData object.
#'
#' @param sce A SingleCellExperiment object to be converted.
#' @param assayName The name of the assay to use as the main data matrix in AnnData. Defaults to "counts".
#' @param useAltExp Logical indicating whether to process and include alternative experiments (altExps). Defaults to TRUE.
#' @return An AnnData object containing the data from the SingleCellExperiment object.
#' @details The function first prints a summary of the input SingleCellExperiment object and checks for alternative experiments (altExps).
#' If altExps are found, they are processed recursively and stored in the `uns` slot of the AnnData object, unless the `useAltExp` parameter is set to FALSE.
#' The specified assay is used as the primary data matrix (anndata.X), and all other assays are stored in the `layers` slot.
#' Dimensional reductions are stored in the `obsm` slot, and metadata is organized into `obs`, `var`, `obsp`, `varp`, and `uns` slots as appropriate.
#' The function includes detailed error handling to ensure that all data is correctly transferred to the AnnData object.
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' library(anndata)
#' sce <- SingleCellExperiment::SingleCellExperiment(list(counts = matrix(1:4, ncol = 2)))
#' ad <- convert_to_anndata(sce)
#' }
#' @importFrom anndata AnnData
#' @export
convert_to_anndata <- function(sce, assayName = "counts", useAltExp = TRUE) {

  # Ensure necessary functions are loaded
  anndata::AnnData

  # Print a summary of the input SCE object
  timestamped_cat("Summary of SingleCellExperiment object:\n\n")
  print(sce)
  cat("\n")
  
  # Process altExperiments
  alt_exps_data <- process_alt_experiments(sce, assayName, useAltExp)
  sce <- alt_exps_data$sce
  alt_exps <- alt_exps_data$alt_exps
  
  # Process the main assay
  X <- process_main_assay(sce, assayName)
  
  # Process other assays
  layers <- process_other_assays(sce, assayName)
  
  # Process dimensional reductions
  obsm <- process_dimensional_reductions(sce)
  
  # Process obs and var data
  obs_data <- extract_data(colData, "obs/colData", sce)
  var_data <- extract_data(rowData, "var/rowData", sce)
  
  # Process metadata and pairwise matrices
  uns <- process_metadata_and_pairwise(sce, alt_exps, X)

  # Create AnnData object
  ad <- AnnData(
    X = X,
    layers = layers,
    obs = obs_data$data,
    var = var_data$data,
    obsm = obsm,
    varm = uns$varm,
    obsp = uns$obsp,
    varp = uns$varp,
    uns = uns$uns
  )

  timestamped_cat("Summary of the AnnData object:\n\n")
  print(ad)
  cat("\n")

  return(ad)
}

#' Convert Seurat or Other Object to SingleCellExperiment
#'
#' This function determines the class of a loaded object and converts it to a SingleCellExperiment object if necessary.
#' It handles Seurat objects, updating old Seurat v2 objects if detected, and converts them to SingleCellExperiment.
#' If the input object is already a SingleCellExperiment, it is returned as is.
#' For other object types, an attempt is made to convert them to SingleCellExperiment.
#'
#' @param data The loaded object to be converted.
#' @return A SingleCellExperiment object.
#' @details The function first checks if the input object is a Seurat object. If it is, it prints a summary of the object
#' and checks for indicators of Seurat v2, updating the object if necessary. It then converts the Seurat object to a SingleCellExperiment.
#' If the object is already a SingleCellExperiment, it is returned directly. For other object types, an attempt is made to convert them to
#' SingleCellExperiment.
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SingleCellExperiment)
#' seurat_obj <- CreateSeuratObject(counts = matrix(1:4, ncol = 2))
#' sce <- convert_seurat_to_sce(seurat_obj)
#' }
#' @import SingleCellExperiment
#' @export
convert_seurat_to_sce <- function(data) {
  object_class <- class(data)
  timestamped_cat("Input object class:", paste(object_class, collapse = ", "), "\n")

  if ("seurat" %in% tolower(object_class)) {
    timestamped_cat("Summary of input Seurat object:\n\n")
    suppressPackageStartupMessages(print(data))
    cat("\n")
    # Use tryCatch to safely check for Seurat v2 or v3 indicators
    raw_data <- tryCatch(
      {
        if (!is.null(data@raw.data)) data@raw.data else NULL
      },
      error = function(e) NULL
    )

    if (!is.null(raw_data)) {
      timestamped_cat("Old Seurat v2 object detected, attempting to update...\n")
      data <- Seurat::UpdateSeuratObject(data)
    }

    # Convert to SingleCellExperiment
    sce <- Seurat::as.SingleCellExperiment(data)
  } else if ("SingleCellExperiment" %in% class(data)) {
    sce <- data
  } else {
    suppressPackageStartupMessages({
      sce <- as(data, "SingleCellExperiment")
    })
  }
  return(sce)
}

#' Command Line Interface for convert2anndata
#'
#' This function serves as a command line interface for the convert2anndata package.
#' It parses command line arguments and calls the appropriate functions to convert
#' a SingleCellExperiment or Seurat object to an AnnData object.
#'
#' @import optparse
#' @import anndata
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
      metavar = "assayName"
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
