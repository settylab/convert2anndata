#' Check and Load Packages
#'
#' This function checks if the required packages are installed, loads them, and suppresses startup messages.
#' If any required packages are not installed, it prompts the user to install them.
#'
#' @param packages A character vector of package names to check and load.
#' @return None. The function either loads the required packages or stops execution if any are missing.
#' @examples
#' check_and_load_packages(c("SingleCellExperiment", "anndata", "optparse"))
check_and_load_packages <- function(packages) {
  missing_packages <- character()
  for (package in packages) {
    suppressPackageStartupMessages({
      if (!require(package, character.only = TRUE, quietly = TRUE)) {
        missing_packages <- c(missing_packages, package)
      }
    })
  }
  if (length(missing_packages) > 0) {
    cat("The following packages are not installed: ", 
        paste(missing_packages, collapse = ", "), "\n")
    cat("Install them by running the following command in R:\n")
    cat("install.packages(c(", paste(sprintf("'%s'", missing_packages), 
        collapse = ", "), "))\n")
    quit(status = 1, save = "no") # Stop execution if any packages are missing
  }
}

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
