#' Extract Data from SingleCellExperiment
#'
#' This function extracts data from the SingleCellExperiment object with careful error handling.
#'
#' @param data_fun A function to extract data (e.g., colData, rowData).
#' @param data_type The type of data being extracted (e.g., "obs/colData", "var/rowData").
#' @param sce A SingleCellExperiment object.
#' @return A list containing the extracted data and internal columns.
#' @importFrom SummarizedExperiment colData rowData
#' @export
extract_data <- function(data_fun, data_type, sce) {
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
