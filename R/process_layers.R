#' Process Layers in Seurat Assay
#'
#' Combines split layers in a Seurat Assay5 object, stores them as non-split layers, and removes the split layers.
#'
#' @param assay_object The Seurat Assay5 object to process.
#' @return A list with split names per layer and updated assay object.
#' @importFrom Seurat GetAssayData SetAssayData
#' @importFrom SeuratObject Layers
#' @importFrom Matrix Matrix
#' @export
process_layers <- function(assay_object) {
  if (!inherits(assay_object, "Assay5")) {
    stop("The assay_object must be an Assay5 object.")
  }

  # Get all layer names
  assay_layers <- Layers(assay_object)
  layer_count <- length(assay_layers)

  if (layer_count == 0) {
    timestamped_cat("No layer found.\n")
    return(list(
      split_names = NULL,
      assay_object = assay_object
    ))
  }

  timestamped_cat("Found", layer_count, "layers:", paste(assay_layers, collapse = ", "), "\n")

  # Parse layer names to get base names and suffixes
  split_info <- strsplit(assay_layers, ".", fixed = TRUE)
  split_groups <- sapply(split_info, `[`, 1) # Base names
  split_map <- split(assay_layers, split_groups)

  # Avoid collision: Exclude non-split layers from split groups
  filtered_split_map <- lapply(names(split_map), function(base_name) {
    group_layers <- split_map[[base_name]]
    split_suffixes <- sapply(strsplit(group_layers, ".", fixed = TRUE), `[`, 2)
    is_split <- !is.na(split_suffixes)
    group_layers[is_split]
  })
  names(filtered_split_map) <- names(split_map) # Preserve names

  # Initialize list to collect split names
  split_names <- list()

  # Process each group
  for (base_name in names(filtered_split_map)) {
    group_layers <- filtered_split_map[[base_name]]
    if (length(group_layers) > 1) {
      # This is a split layer, concatenate the parts
      timestamped_cat("Combining split layers for base name:", base_name, "\n")
      layer_counts_list <- lapply(group_layers, function(layer_name) {
        layer_counts <- GetAssayData(assay_object, layer = layer_name)
        layer_counts <- layer_counts[rownames(assay_object), , drop = FALSE] # Align rows
        return(layer_counts)
      })
      combined_layer_counts <- do.call(cbind, layer_counts_list)

      # Convert to CsparseMatrix format
      combined_layer_counts <- methods::as(combined_layer_counts, "CsparseMatrix")

      # Extract split names
      group_names <- sapply(strsplit(group_layers, ".", fixed = TRUE), `[`, 2)
      split_names[[base_name]] <- rep(group_names, sapply(layer_counts_list, ncol))

      # Order columns to match assay
      original_order <- colnames(assay_object)
      combined_columns <- colnames(combined_layer_counts)
      column_order <- match(original_order, combined_columns)
      combined_layer_counts <- combined_layer_counts[, column_order, drop = FALSE]
      split_names[[base_name]] <- split_names[[base_name]][order(column_order)]

      # Determine the name for the combined layer
      new_layer_name <- base_name
      if (new_layer_name %in% assay_layers) {
        suffix <- 1
        while (paste0(new_layer_name, ".", suffix) %in% assay_layers) {
          suffix <- suffix + 1
        }
        new_layer_name <- paste0(new_layer_name, ".", suffix)
      }

      # Store combined layer in the assay object
      assay_object <- tryCatch(
        {
          SetAssayData(assay_object, layer = new_layer_name, new.data = combined_layer_counts)
        },
        error = function(e) {
          stop("Failed to set combined layer for base name:", base_name, "\n", e$message)
        }
      )

      # Remove the original split layers
      for (layer_name in group_layers) {
        tryCatch(
          {
            assay_object <- SetAssayData(assay_object, layer = layer_name, new.data = NULL)
          },
          error = function(e) {
            timestamped_cat("Error removing split layer:", layer_name, "\n")
          }
        )
      }
    }
  }

  list(
    split_names = split_names,
    assay_object = assay_object
  )
}
