#' Extract Pairwise Data from SingleCellExperiment
#'
#' This function extracts pairwise data (colPairs or rowPairs) from a SingleCellExperiment
#' object and returns it as a list of sparse matrices.
#'
#' @param pairs_function A function to extract pairwise data (e.g., colPairs or rowPairs).
#' @param sce A SingleCellExperiment object from which to extract the pairwise data.
#' @return A list of pairwise data as sparse matrices.
#' @importFrom SingleCellExperiment colPairs rowPairs
#' @importFrom Matrix sparseMatrix
#' @importFrom methods is
#' @export
extract_pairs <- function(pairs_function, sce) {
  pairs_list <- list()
  pairs <- pairs_function(sce)

  if (length(pairs) == 0) {
    timestamped_cat("No pairwise data found.\n")
    return(pairs_list)
  }

  if ("index" %in% names(pairs) && "distance" %in% names(pairs) &&
      is(pairs[["index"]], "SelfHits") && is(pairs[["distance"]], "SelfHits")) {
    timestamped_cat("Detected index and distance pairwise data in SelfHit format.\n")
    index_hits <- pairs[["index"]]
    distance_hits <- pairs[["distance"]]

    if (all(queryHits(index_hits) == queryHits(distance_hits)) && all(subjectHits(index_hits) == subjectHits(distance_hits))) {
      timestamped_cat(
        "Index and distance 'from' and 'to' columns agree. Combining info to",
        "create 'distance' and 'connectivity' sparse natrix.\n"
      )
      pairs_list[["distance"]] <- sparseMatrix(
        i = queryHits(index_hits),
        j = mcols(index_hits)[[1]],
        x = mcols(distance_hits)[[1]],
        dims = c(ncol(sce), ncol(sce))
      )
      pairs_list[["connectivity"]] <- sparseMatrix(
        i = queryHits(index_hits),
        j = mcols(index_hits)[[1]],
        x = rep(1, length(index_hits)),
        dims = c(ncol(sce), ncol(sce))
      )
      pairs[["index"]] <- NULL
      pairs[["distance"]] <- NULL
    } else {
      timestamped_cat(
        "WARNING: Index and distance 'from' and 'to' columns do not agree.",
        "Using them individually and intereting 'to' column as cell index.\n"
      )
    }
  }

  for (name in names(pairs)) {
    if (is(pairs[[name]], "SelfHits")) {
      if (length(mcols(pairs[[name]])) == 1) {
        pairs_list[[name]] <- sparseMatrix(
          i = queryHits(pairs[[name]]),
          j = subjectHits(pairs[[name]]),
          x = mcols(pairs[[name]])[[1]],
          dims = c(ncol(sce), ncol(sce))
        )
        timestamped_cat("Created sparse matrix for", name, "paired data.\n")
      } else if (length(mcols(pairs[[name]])) > 1) {
        for (mcol_name in colnames(mcols(pairs[[name]]))) {
          pairs_list[[paste0(name, "_", mcol_name)]] <- sparseMatrix(
            i = queryHits(pairs[[name]]),
            j = subjectHits(pairs[[name]]),
            x = mcols(pairs[[name]])[[mcol_name]],
            dims = c(ncol(sce), ncol(sce))
          )
          timestamped_cat("Created sparse matrix for", name, "_", mcol_name, "paired data.\n")
        }
      } else {
        pairs_list[[name]] <- sparseMatrix(
          i = queryHits(pairs[[name]]),
          j = subjectHits(pairs[[name]]),
          x = rep(1, length(pairs[[name]])),
          dims = c(ncol(sce), ncol(sce))
        )
        timestamped_cat("Created sparse matrix for", name, "paired data with default values.\n")
      }
    } else {
      pairs_list[[name]] <- as(pairs[[name]], "dgCMatrix")
      timestamped_cat("Converted", name, "paired data to dgCMatrix.\n")
    }
  }
  

  return(pairs_list)
}
