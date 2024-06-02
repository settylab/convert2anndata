#' Process Alternative Experiments
#'
#' This function processes alternative experiments (altExps) in the SingleCellExperiment object.
#' If altExps are found, they are processed recursively and stored in the `uns` slot of the AnnData object.
#'
#' @param sce A SingleCellExperiment object.
#' @param assayName The name of the assay to use as the main data matrix in AnnData.
#' @param useAltExp Logical indicating whether to process and include alternative experiments (altExps).
#' @return A list containing the updated SingleCellExperiment object and the altExps data.
#' @importFrom SingleCellExperiment altExpNames altExp removeAltExps
#' @export
process_alt_experiments <- function(sce, assayName, useAltExp) {
  timestamped_cat("Checking for altExperiments...\n")
  alt_exps <- altExpNames(sce)
  alt_exps_data <- list()

  if (length(alt_exps) > 0) {
    timestamped_cat(
      "WARNING: The SingleCellExperiment object contains",
      "alternative experiments (altExps) which cannot ideally be",
      "reflected in the AnnData object.\n"
    )
    if (isTRUE(useAltExp)) {
      timestamped_cat(
        "The following altExps will be processed and stored in",
        "uns['altExperiments']: ",
        paste(alt_exps, collapse = ", "), "\n"
      )
      for (alt_exp in alt_exps) {
        timestamped_cat(sprintf("Starting processing of altExp '%s'...\n", alt_exp))
        alt_sce <- altExp(sce, alt_exp)
        alt_assayName <- assayName
        if (!(assayName %in% names(assays(alt_sce)))) {
          alt_assayName <- names(assays(alt_sce))[1]
          timestamped_cat(sprintf(
            "WARNING: The specified assay '%s' is not ",
            "available in altExp '%s'. Using '%s' instead.\n",
            assayName, alt_exp, alt_assayName
          ))
        }
        alt_ad <- convert_to_anndata(alt_sce, alt_assayName)
        alt_exps_data[[alt_exp]] <- alt_ad
        timestamped_cat(sprintf("Processed altExp '%s'.\n", alt_exp))
      }
    } else {
      timestamped_cat(
        "Recursive recovery of altExperiments is disabled,",
        "removing experiments:",
        paste(alt_exps, collapse = ", "), "\n"
      )
    }
    sce <- removeAltExps(sce)
  }

  return(list(sce = sce, alt_exps = alt_exps_data))
}
