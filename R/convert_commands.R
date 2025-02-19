#' Convert SeuratCommand objects into a serializable format
#'
#' This function takes a list of SeuratCommand objects and converts them into
#' JSON-serializable lists, preserving command name, parameters, timestamps,
#' and the assay used.
#'
#' @param commands_list A list of SeuratCommand objects from Seurat metadata.
#' @return A list of simplified, serializable command information.
#' @export
convert_commands <- function(commands_list) {
  lapply(commands_list, function(cmd) {
    params_serializable <- lapply(cmd@params, function(param) {
      if (is.function(param)) {
        return(deparse(param)[1])  # Store function as a string
      } else {
        return(param)
      }
    })

    list(
      command = cmd@call.string,  # Extract function call as string
      time = format(cmd@time.stamp, "%Y-%m-%d %H:%M:%S"),  # Convert timestamp
      assay_used = cmd@assay.used,  # Include the assay name
      params = params_serializable  # Filter out non-serializable params
    )
  })
}