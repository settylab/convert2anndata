#' Convert Seurat or Other Object to SingleCellExperiment
#'
#' This function converts a loaded object to a `SingleCellExperiment` object if necessary.
#' It first attempts to use Seurat's built-in conversion function. If this fails (e.g., due to multiple layers),
#' it performs a custom conversion, preserving multiple assays, paired data (such as distance matrices),
#' and handling mismatches appropriately. It also attempts to transfer unstructured data by storing it in
#' the `metadata` slot of the SingleCellExperiment object. If the input object is already a `SingleCellExperiment`,
#' it is returned as is.
#'
#' For assays with multiple layers (e.g., `Assay5` in Seurat v5), the function processes split layers,
#' combines them where necessary, and records their names in the metadata. Layers with mismatched sizes or
#' dimensions are moved to `altExps`, with name collisions resolved by appending unique suffixes.
#'
#' Dimensional reductions, graphs, feature metadata, and unstructured data are also transferred to ensure
#' compatibility with `SingleCellExperiment`.
#'
#' @param data The loaded object to be converted.
#' @return A `SingleCellExperiment` object.
#' @importFrom Seurat as.SingleCellExperiment DefaultAssay
#' @importFrom SingleCellExperiment SingleCellExperiment altExp<- altExpNames colPair<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom methods as slotNames slot
#' @export
convert_seurat_to_sce <- function(data) {
  object_class <- class(data)
  timestamped_cat("Input object class:", paste(object_class, collapse = ", "), "\n")

  if ("seurat" %in% tolower(object_class)) {
    timestamped_cat("Converting Seurat object...\n")

    # Update Seurat object if necessary
    data <- update_seurat_object(data)

    # Attempt to use Seurat's built-in conversion function (disabled as requested)
    if (FALSE) {
      timestamped_cat("Attempting to convert using Seurat's built-in function...\n")
      sce <- tryCatch(
        {
          sce_obj <- as.SingleCellExperiment(data)
          timestamped_cat("Conversion using Seurat's built-in function succeeded.\n")
          return(sce_obj)
        },
        error = function(e) {
          timestamped_cat("Seurat's built-in conversion failed: ", e$message, "\n")
          NULL
        }
      )
      if (!is.null(sce)) {
        return(sce)
      }
    }

    # Proceeding with custom conversion
    timestamped_cat("Proceeding with custom conversion...\n")
    if (nrow(data) == 0 || ncol(data) == 0) {
      stop("The Seurat object contains no data. Ensure the object has non-zero features and cells.")
    }

    assays_list <- list()
    altExp_list <- list()
    assays_to_move <- character(0)

    # Extract all assay names
    assay_names <- names(data@assays)
    timestamped_cat("Found assays:", paste(assay_names, collapse = ", "), "\n")

    # Extract reference feature and cell names from the default assay
    default_assay_name <- DefaultAssay(data)
    default_assay <- data@assays[[default_assay_name]]
    ref_cells <- colnames(default_assay)
    ref_features <- rownames(default_assay)

    # Collect metadata
    metadata_data <- data@meta.data

    # Process each assay
    for (assay_name in assay_names) {
      assay_object <- data@assays[[assay_name]]
      timestamped_cat("Processing assay:", assay_name, "\n")

      if (inherits(assay_object, "Assay5")) {
        # Handle multi-layered assays
        layer_result <- process_layers(assay_object)
        split_names <- layer_result$split_names

        # Store split layer information in metadata
        if (!is.null(split_names)) {
          for (base_name in names(split_names)) {
            metadata_data[[paste0(assay_name, "_", base_name, "_split")]] <- split_names[[base_name]]
          }
        }

        # Add layers to assays list or move mismatched layers to altExp
        for (layer_name in Layers(layer_result$assay_object)) {
          layer_data <- GetAssayData(layer_result$assay_object, layer = layer_name)
          if (!all(rownames(layer_data) == ref_features) || !all(colnames(layer_data) == ref_cells)) {
            alt_name <- paste0(assay_name, "_", layer_name)
            while (alt_name %in% names(altExp_list)) {
              alt_name <- paste0(alt_name, "_alt")
            }
            altExp_list[[alt_name]] <- SingleCellExperiment(assays = list(counts = layer_data))
            timestamped_cat("Moved layer", layer_name, "to altExp due to mismatched cells/features:", alt_name, "\n")
          } else {
            assays_list[[paste0(assay_name, "_", layer_name)]] <- layer_data
          }
        }
      } else {
        # Process standard assays
        counts_matrix <- extract_counts_matrix(assay_object)
        if (!all(rownames(counts_matrix) == ref_features) || !all(colnames(counts_matrix) == ref_cells)) {
          # Move mismatched assay to altExp
          alt_name <- assay_name
          while (alt_name %in% names(altExp_list)) {
            alt_name <- paste0(alt_name, "_alt")
          }
          altExp_list[[alt_name]] <- SingleCellExperiment(assays = list(counts = counts_matrix))
          timestamped_cat("Moved assay", assay_name, "to altExp due to mismatched cells/features:", alt_name, "\n")
        } else {
          assays_list[[assay_name]] <- counts_matrix
        }
      }
    }

    # Create the SingleCellExperiment object
    sce <- SingleCellExperiment(assays = assays_list, colData = metadata_data)

    # Add altExps
    for (alt_name in names(altExp_list)) {
      altExp(sce, alt_name) <- altExp_list[[alt_name]]
    }

    # Transfer dimensional reductions
    sce <- preserve_reductions_seurat(data, sce)

    # Transfer graphs (e.g., distance matrices)
    if (length(data@graphs) > 0) {
      timestamped_cat("Found", length(data@graphs), "graphs. Moving to colPair.")
      for (graph_name in names(data@graphs)) {
        graph <- data@graphs[[graph_name]]
        if (all(rownames(graph) %in% colnames(sce)) && all(colnames(graph) %in% colnames(sce))) {
          colPair(sce, graph_name) <- convert_graph_to_colPair(graph, colnames(sce))
        } else {
          metadata(sce)[[paste0("graph_", graph_name)]] <- graph
        }
      }
    }

    # Transfer miscellaneous and unstructured data
    known_slots <- c(
      "assays", "meta.data", "active.assay", "active.ident", "graphs", "reductions",
      "project.name" # , "misc", "commands", "version", "tools"
    )
    additional_slots <- setdiff(slotNames(data), known_slots)
    for (slot_name in additional_slots) {
      slot_data <- slot(data, slot_name)
      if (!is.null(slot_data) && length(slot_data) > 0) {
        metadata(sce)[[slot_name]] <- slot_data
      }
    }

    timestamped_cat("Custom conversion to SingleCellExperiment completed successfully.\n")
    return(sce)
  } else if ("SingleCellExperiment" %in% object_class) {
    timestamped_cat("Input object is already a SingleCellExperiment. Returning as is.\n")
    return(data)
  } else {
    # Attempt to coerce other object types to SingleCellExperiment
    suppressPackageStartupMessages({
      sce <- tryCatch(
        {
          as(data, "SingleCellExperiment")
        },
        error = function(e) {
          timestamped_cat("Error converting to SingleCellExperiment: ", e$message, "\n")
          stop("Failed to convert object to SingleCellExperiment: ", e$message)
        }
      )
    })
    timestamped_cat("Conversion to SingleCellExperiment completed successfully.\n")
    return(sce)
  }
}
