#' Command Line Interface for convert2anndata
#'
#' Bidirectional CLI: dispatches on the input file extension. An `.rds` input
#' is treated as a Seurat or SingleCellExperiment object and converted to an
#' AnnData `.h5ad` file. An `.h5ad` input is read with `anndata::read_h5ad()`
#' and converted to a Seurat (default) or SingleCellExperiment object saved
#' to `.rds`.
#'
#' @import optparse
#' @importFrom anndata write_h5ad read_h5ad
#' @importFrom optparse make_option OptionParser parse_args
#' @export
cli_convert <- function() {
  anndata::write_h5ad
  anndata::read_h5ad

  description <- paste(
    "Bidirectional converter between AnnData (.h5ad) and Seurat / SingleCellExperiment (.rds).",
    "Input format is detected by extension: .rds -> .h5ad, .h5ad -> .rds."
  )

  option_list <- list(
    make_option(c("-i", "--input"),
      type = "character", default = NULL,
      help = "Path to the input file (.rds or .h5ad). Required.",
      metavar = "file"
    ),
    make_option(c("-o", "--output"),
      type = "character", default = NULL,
      help = paste(
        "Path to the output file. If not specified, the output path is derived",
        "by swapping the extension (.rds <-> .h5ad)."
      ),
      metavar = "file"
    ),
    make_option(c("-a", "--assay"),
      type = "character", default = "counts",
      help = paste(
        "For .rds -> .h5ad: assay to use as anndata.X (defaults to 'counts').",
        "For .h5ad -> .rds: layer name to use as the counts assay (defaults to 'counts')."
      ),
      metavar = "assay_name"
    ),
    make_option(c("-d", "--disable-recursive-altExp"),
      action = "store_true", default = FALSE,
      help = "Disable recursive recovery of altExperiments and discard them instead.",
      metavar = "boolean"
    ),
    make_option(c("-t", "--target"),
      type = "character", default = "seurat",
      help = paste(
        "For .h5ad -> .rds only: target object type, 'seurat' (default) or 'sce'."
      ),
      metavar = "type"
    )
  )

  opt_parser <- OptionParser(option_list = option_list, description = description)
  opt <- parse_args(opt_parser)

  if (is.null(opt$input)) {
    stop("No input file provided. Use --input to specify the input file.", call. = FALSE)
  }

  is_h5ad <- grepl("\\.h5ad$", opt$input, ignore.case = TRUE)
  is_rds  <- grepl("\\.[rR][dD][sS]$", opt$input)

  if (!is_h5ad && !is_rds) {
    stop("Input file must end in .rds or .h5ad. Got: ", opt$input, call. = FALSE)
  }

  if (is_rds) {
    if (is.null(opt$output)) {
      opt$output <- sub("\\.[rR][dD][sS]$", ".h5ad", opt$input, ignore.case = TRUE)
    }
    timestamped_cat("Loading data from:", opt$input, "\n")
    data <- readRDS(opt$input)
    sce <- convert_seurat_to_sce(data)
    timestamped_cat("Data loaded and converted successfully if needed.\n")
    ad <- convert_to_anndata(sce, opt$assay, useAltExp = !isTRUE(opt$`disable-recursive-altExp`))
    timestamped_cat("Saving the AnnData object to:", opt$output, "\n")
    write_h5ad(ad, opt$output)
    timestamped_cat("Conversion complete:", opt$output, "\n")
  } else {
    if (is.null(opt$output)) {
      opt$output <- sub("\\.h5ad$", ".rds", opt$input, ignore.case = TRUE)
    }
    target <- tolower(opt$target)
    if (!target %in% c("seurat", "sce")) {
      stop("--target must be 'seurat' or 'sce'. Got: ", opt$target, call. = FALSE)
    }
    timestamped_cat("Loading AnnData from:", opt$input, "\n")
    adata <- anndata::read_h5ad(opt$input)
    if (target == "seurat") {
      out <- convert_anndata_to_seurat(adata, counts_layer = opt$assay)
    } else {
      out <- convert_anndata_to_sce(adata, counts_layer = opt$assay)
    }
    timestamped_cat("Saving", target, "object to:", opt$output, "\n")
    saveRDS(out, opt$output)
    timestamped_cat("Conversion complete:", opt$output, "\n")
  }
}
