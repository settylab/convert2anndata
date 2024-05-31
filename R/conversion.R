#' Print Messages with Timestamp
#'
#' This function prints messages with a timestamp.
#'
#' @param ... The messages to print.
#' @return None. The function prints the messages to the console.
#' @examples
#' timestamped_cat("This is a message.")
timestamped_cat <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ...)
}

#' Convert SingleCellExperiment to AnnData
#'
#' This function converts a SingleCellExperiment (SCE) object to an AnnData object.
#' It processes assays, dimensional reductions, and metadata, including alternative experiments (altExps).
#' The main assay is used as the primary data matrix in the AnnData object.
#'
#' @param sce A SingleCellExperiment object to be converted.
#' @param assay_name The name of the assay to use as the main data matrix in AnnData. Defaults to "counts".
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
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import anndata
#' @export
convert_to_anndata <- function(sce, assay_name = "counts", useAltExp = TRUE) {

  # Ensure SingleCellExperiment functions are loaded
  SingleCellExperiment::altExpNames
  SingleCellExperiment::altExp
  SingleCellExperiment::removeAltExps
  SummarizedExperiment::assay
  S4Vectors::metadata
  anndata::AnnData

  # Print a summary of the input SCE object
  timestamped_cat("Summary of SingleCellExperiment object:\n\n")
  print(sce)
  cat("\n")
  
  timestamped_cat("Checking for altExperiments...\n")
  alt_exps <- altExpNames(sce)
  alt_exps_data = list()
  if (length(alt_exps) > 0) {
    timestamped_cat("WARNING: The SingleCellExperiment object contains",
                    "alternative experiments (altExps) which cannot ideally be",
                    "reflected in the AnnData object.\n")
    if (isTRUE(useAltExp)) {
      timestamped_cat("The following altExps will be processed and stored in",
                      "uns['altExperiments']: ",
                      paste(alt_exps, collapse = ", "), "\n")
      for (alt_exp in alt_exps) {
        timestamped_cat(sprintf("Starting processing of altExp '%s'...\n", alt_exp))
        alt_sce <- altExp(sce, alt_exp)
        alt_assay_name <- assay_name
        if (!(assay_name %in% names(assays(alt_sce)))) {
          alt_assay_name <- names(assays(alt_sce))[1]
          timestamped_cat(sprintf("WARNING: The specified assay '%s' is not ",
                                  "available in altExp '%s'. Using '%s' instead.\n",
                                  assay_name, alt_exp, alt_assay_name))
        }
        alt_ad <- convert_to_anndata(alt_sce, alt_assay_name)
        alt_exps_data[[alt_exp]] <- alt_ad
        timestamped_cat(sprintf("Processed altExp '%s'.\n", alt_exp))
      }
    } else {
      timestamped_cat("Recursive recovery of altExperiments is disabled,",
                      "removing experiments:",
                      paste(alt_exps, collapse = ", "), "\n")
    }
    sce <- removeAltExps(sce)
  }

  timestamped_cat(sprintf("Processing assay '%s' for anndata.X...\n", 
                          assay_name))
  if (!(assay_name %in% names(assays(sce)))) {
    timestamped_cat(
      sprintf(
        "Error: The specified assay '%s' is not available in the provided ",
        "SingleCellExperiment object.",
        assay_name
      ),
      "Use -a or --assay to specify an available assay.\n"
    )
    timestamped_cat("Available assays are:", 
                    paste(names(assays(sce)), collapse = ", "), "\n")
    quit(status = 1, save = "no")
  }
  X <- Matrix::t(assay(sce, assay_name))
  timestamped_cat(sprintf("Using '%s' assay as the main data matrix.\n", 
                          assay_name))
  
  timestamped_cat(sprintf("Processing assays other than '%s'...\n", assay_name))
  all_assays <- assays(sce)
  all_assays <- all_assays[!names(all_assays) %in% assay_name]
  all_assays <- lapply(all_assays, Matrix::t)
  timestamped_cat("Assays processed.\n")

  timestamped_cat("Processing dimensional reductions...\n")
  available_reductions <- names(reducedDims(sce))
  if (length(available_reductions) > 0) {
    obsm <- lapply(available_reductions, function(rd_name) {
      reducedDim(sce, rd_name)
    })
    names(obsm) <- paste("X", tolower(available_reductions), sep = "_")
    timestamped_cat("Dimensional reductions processed.\n")
  } else {
    obsm <- list()
    timestamped_cat("No dimensional reductions found.\n")
  }
  
  timestamped_cat("Processing and filtering the obs/colData and var/rowData...\n")
  
  extract_data <- function(data_fun, data_type) {
    internal_columns <- NULL
    tryCatch({
      internal_columns <- colnames(data_fun(sce, internal = TRUE))
      data <- as.data.frame(data_fun(sce, internal = TRUE))
      timestamped_cat(sprintf("Successfully extracted %s with internal=TRUE.\n", 
                              data_type))
      return(list(data = data, internal_columns = internal_columns))
    }, error = function(e) {
      timestamped_cat(sprintf("WARNING: Failed to extract %s with internal=TRUE.\n", 
                              data_type))
      tryCatch({
        data <- as.data.frame(data_fun(sce))
        timestamped_cat(sprintf("Successfully extracted %s without internal=TRUE.\n", 
                                data_type))
        missing_columns <- setdiff(internal_columns, colnames(data))
        if (length(missing_columns) > 0) {
          timestamped_cat("WARNING: The following columns are missing without ",
                          "internal=TRUE and some meta information might be lost: ",
                          paste(missing_columns, collapse = ", "), "\n")
        }
        return(list(data = data, internal_columns = internal_columns))
      }, error = function(e) {
        timestamped_cat(sprintf("ERROR: Failed to extract %s.\n", data_type))
        stop(sprintf("Unable to extract %s from the SingleCellExperiment object.", 
                     data_type), call. = FALSE)
      })
    })
  }
  
  # Extract obs_data and var_data with careful error handling
  obs_result <- extract_data(colData, "obs/colData")
  obs_data <- obs_result$data
  
  var_result <- extract_data(rowData, "var/rowData")
  var_data <- var_result$data
  
  reduction_prefixes <- paste0(tolower(available_reductions), "\\.")
  reduction_columns <- grep("^reducedDims\\.", names(obs_data), value = TRUE)
  if (length(reduction_columns) > 0) {
    obs_data <- obs_data[, !(names(obs_data) %in% reduction_columns), 
                         drop = FALSE]
  }
  timestamped_cat("obs/colData processed.\n")
  
  timestamped_cat("Gathering metadata and pairwise matrices...\n")
  uns <- c(metadata(sce), int_metadata(sce))
  obsp <- list()
  varp <- list()
  varm <- list()
  n_obs <- nrow(X)
  n_var <- ncol(X)
  for (name in names(uns)) {
    item <- uns[[name]]
    if (is.matrix(item) && all(dim(item) == c(n_obs, n_obs))) {
      obsp[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to obsp
      timestamped_cat(sprintf("Transferred '%s' to pairwise observations (obsp).\n", 
                              name))
    } else if (is.matrix(item) && nrow(item) == n_obs) {
      obsm[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to obsm
      timestamped_cat(sprintf("Transferred '%s' to observation matrices (obsm).\n", 
                              name))
    }
    if (is.matrix(item) && all(dim(item) == c(n_var, n_var))) {
      varp[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to obsp
      timestamped_cat(sprintf("Transferred '%s' to pairwise variables (varp).\n", 
                              name))
    } else if (is.matrix(item) && nrow(item) == n_var) {
      varm[[name]] <- item
      uns[[name]] <- NULL # Remove from uns if added to varm
      timestamped_cat(sprintf("Transferred '%s' to variable matrices (varm).\n", 
                              name))
    }
  }
  if (length(alt_exps_data) > 0) {
    uns[["altExperiments"]] <- alt_exps_data
  }
  timestamped_cat("Metadata and pairwise data organized.\n")
  
  timestamped_cat("Making AnnData object...\n")
  ad <- AnnData(
    X = X,
    layers = all_assays,
    obs = obs_data,
    var = var_data,
    obsm = obsm,
    varm = varm,
    obsp = obsp,
    varp = varp,
    uns = uns
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
  if ("seurat" %in% tolower(class(data))) {
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
