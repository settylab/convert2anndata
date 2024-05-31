#' Print Messages with Timestamp
#'
#' This function prints messages with a timestamp.
#'
#' @param ... The messages to print.
#' @return None. The function prints the messages to the console.
#' @examples
#' timestamped_cat("This is a message.")
#' @export
timestamped_cat <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ...)
}
