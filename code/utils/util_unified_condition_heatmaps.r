# util_unified_condition_heatmaps.r - ARCHIVED VERSION
# Creates a single heatmap per timepoint containing all patients.
# This version is modified for the self-contained summer report archive.
# Key change: All file paths now point to a local `./data` directory.

here::i_am("code/utils/util_unified_condition_heatmaps.r")

library(here)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(readr)

#' Load sample metadata with GMM information
#' @param dataset_name Name of the dataset
#' @param run_name Run identifier
#' @return Data frame with sample metadata
load_sample_metadata_unified <- function(dataset_name = "HIV_PBMC", run_name = "6-20k") {

  # MODIFIED PATH: Point to the local data directory
  data_dir <- here::here("data", "gloscope_results", run_name)

  # Look for GMM parameter files
  param_files <- list.files(data_dir, pattern = paste0("gmm_params_.*_defaultGMM_", run_name, "\\.rds"), full.names = TRUE)

  if (length(param_files) == 0) {
    stop("No GMM parameter files found in ", data_dir)
  }

  # Extract sample information
  sample_ids <- gsub(paste0(".*gmm_params_(.*)_defaultGMM_", run_name, "\\.rds"), "\\1", basename(param_files))

  sample_info <- data.frame(
    sample_id = sample_ids,
    stringsAsFactors = FALSE
  )

  # Parse patient and timepoint
  sample_info$patient <- gsub("_.*", "", sample_info$sample_id)
  sample_info$timepoint <- gsub(".*_", "", sample_info$sample_id)

  # Standardize timepoint names
  sample_info$timepoint_clean <- sample_info$timepoint
  sample_info$timepoint_clean <- gsub("Weeks?", " Weeks", sample_info$timepoint_clean)
  sample_info$timepoint_clean <- gsub("Months?", " Months", sample_info$timepoint_clean)
  sample_info$timepoint_clean <- gsub("Years?", " Year", sample_info$timepoint_clean)
  sample_info$timepoint_clean <- gsub("PreInfection", "Pre-Infection", sample_info$timepoint_clean)
  # Trim whitespace to prevent issues like "6  Months"
  sample_info$timepoint_clean <- trimws(gsub("\\s+", " ", sample_info$timepoint_clean))

  # Load GMM parameters and K values
  sample_info$k_components <- NA
  sample_info$model_type <- NA

  for (i in 1:nrow(sample_info)) {
    # MODIFICATION: Use data_dir
    param_file <- file.path(data_dir, paste0("gmm_params_", sample_info$sample_id[i], "_defaultGMM_", run_name, ".rds"))

    if (file.exists(param_file)) {
      tryCatch({
        params <- readRDS(param_file)
        sample_info$k_components[i] <- ncol(params$mu)
        sample_info$model_type[i] <- params$modelName %||% "Unknown"
      }, error = function(e) {
        cat("Warning: Could not load parameters for", sample_info$sample_id[i], "\n")
      })
    }
  }

  sample_info <- sample_info %>%
    arrange(timepoint_clean, patient)

  return(sample_info)
}

#' Load cell-type anchors from seurat object
#' @param dataset_name Name of the dataset
#' @return Data frame with cell-type anchors in PC space
load_celltype_anchors_unified <- function(dataset_name = "HIV_PBMC") {

  # MODIFIED PATH: Point to the local data directory
  seurat_file <- here::here("data", "processed_data", "default_seurat.rds")

  if (!file.exists(seurat_file)) {
    warning("Seurat file not found: ", seurat_file)
    return(NULL)
  }

  seurat_obj <- readRDS(seurat_file)

  # Extract PCA coordinates and metadata
  pca_coords <- seurat_obj[["pca"]]@cell.embeddings[, 1:10]
  colnames(pca_coords) <- paste0("PC", 1:10)

  metadata <- seurat_obj@meta.data

  pc_data <- as.data.frame(pca_coords) %>%
    bind_cols(metadata)

  # Calculate anchors (mean PC coordinates per cell type)
  anchors <- pc_data %>%
    group_by(cell_type) %>%
    summarise(across(PC1:PC10, ~ mean(.x, na.rm = TRUE)), .groups = 'drop')

  return(anchors)
}

#' Create unified mu matrix for a single timepoint with anchors
#' @param dataset_name Name of the dataset
#' @param timepoint_name Name of the timepoint
#' @param sample_metadata Data frame with sample metadata
#' @param run_name Run identifier
#' @param use_pcs Number of PC dimensions to use
#' @return List containing mu matrix and sample information
create_unified_mu_matrix <- function(dataset_name, timepoint_name, sample_metadata,
                                   run_name = "6-20k", use_pcs = 10) {

  # Get samples for this timepoint
  timepoint_samples <- sample_metadata %>%
    filter(timepoint_clean == timepoint_name)
    
  if (nrow(timepoint_samples) == 0) {
    warning("No samples found for timepoint: ", timepoint_name)
    return(NULL)
  }

  cat("Processing timepoint:", timepoint_name, "with", nrow(timepoint_samples), "samples\n")

  # Initialize containers
  mu_list <- list()
  sample_info_list <- list()

  # Load mu vectors for each sample
  for (i in 1:nrow(timepoint_samples)) {
    sample_id <- timepoint_samples$sample_id[i]
    patient <- timepoint_samples$patient[i]
    k_comp <- timepoint_samples$k_components[i]

    cat("  Loading", sample_id, "(", patient, ", K =", k_comp, ")\n")

    # Load GMM parameters
    # MODIFIED PATH: Point to the local data directory
    param_file <- here::here("data", "gloscope_results", run_name,
                             paste0("gmm_params_", sample_id, "_defaultGMM_", run_name, ".rds"))

    if (file.exists(param_file)) {
      tryCatch({
        params <- readRDS(param_file)

        if (!is.null(params$mu)) {
          mu_matrix <- params$mu

          # Use specified number of PCs
          if (nrow(mu_matrix) > use_pcs) {
            mu_matrix <- mu_matrix[1:use_pcs, ]
          }

          # Add each component to the list
          for (comp in 1:ncol(mu_matrix)) {
            unique_id <- paste0(sample_id, "_C", comp)
            mu_list[[unique_id]] <- mu_matrix[, comp]

            sample_info_list[[unique_id]] <- data.frame(
              unique_id = unique_id,
              sample_id = sample_id,
              component = comp,
              patient = patient,
              timepoint = timepoint_name,
              is_anchor = FALSE,
              stringsAsFactors = FALSE
            )
          }
        }
      }, error = function(e) {
        cat("  Warning: Could not load mu vectors for", sample_id, "\n")
      })
    }
  }

  if (length(mu_list) == 0) {
    cat("No valid mu vectors found for", timepoint_name, "\n")
    return(NULL)
  }

  # Load cell-type anchors
  anchors <- load_celltype_anchors_unified(dataset_name)

  if (!is.null(anchors)) {
    # Add anchors to mu_list
    for (a in 1:nrow(anchors)) {
      celltype <- anchors$cell_type[a]
      anchor_id <- paste0("CellType-", celltype)
      mu_list[[anchor_id]] <- as.numeric(anchors[a, paste0("PC", 1:use_pcs)])

      sample_info_list[[anchor_id]] <- data.frame(
        unique_id = anchor_id,
        sample_id = anchor_id,
        component = 0,
        patient = "Anchor",
        timepoint = timepoint_name,
        is_anchor = TRUE,
        stringsAsFactors = FALSE
      )
    }
    cat("  Added", nrow(anchors), "cell-type anchors\n")
  }

  # Convert to matrix
  mu_matrix <- do.call(rbind, mu_list)
  sample_info_df <- do.call(rbind, sample_info_list)
  rownames(mu_matrix) <- sample_info_df$unique_id

  cat("  Final matrix dimensions:", nrow(mu_matrix), "entries\n")

  return(list(
    mu_matrix = mu_matrix,
    sample_info = sample_info_df
  ))
}

#' Create unified heatmap for a single timepoint
#' @param dataset_name Name of the dataset
#' @param timepoint_name Name of the timepoint
#' @param sample_metadata Data frame with sample metadata
#' @param run_name Run identifier
#' @param use_pcs Number of PC dimensions to use
#' @return File path of generated heatmap
create_unified_timepoint_heatmap <- function(dataset_name, timepoint_name, sample_metadata,
                                           run_name = "6-20k", use_pcs = 10) {

  # Create unified mu matrix with anchors
  mu_data <- create_unified_mu_matrix(dataset_name, timepoint_name, sample_metadata, run_name, use_pcs)

  if (is.null(mu_data)) {
    warning("Could not create mu matrix for timepoint: ", timepoint_name)
    return(NULL)
  }

  # Calculate Euclidean distance matrix
  dist_matrix <- as.matrix(dist(mu_data$mu_matrix, method = "euclidean"))

  # Create annotation dataframe
  sample_info <- mu_data$sample_info

  annotation_df <- data.frame(
    Sample = sample_info$sample_id,
    Type = ifelse(sample_info$is_anchor, "Anchor", "μ-vector"),
    row.names = sample_info$unique_id
  )

  # Create color scheme
  unique_samples <- unique(sample_info$sample_id[!sample_info$is_anchor])

  if (length(unique_samples) <= 8) {
    sample_colors <- RColorBrewer::brewer.pal(max(3, length(unique_samples)), "Set2")
  } else {
    sample_colors <- rainbow(length(unique_samples))
  }
  names(sample_colors) <- unique_samples

  # Add anchor colors (black for all anchors)
  anchor_ids <- unique(sample_info$sample_id[sample_info$is_anchor])
  if (length(anchor_ids) > 0) {
    anchor_colors <- rep("#000000", length(anchor_ids))
    names(anchor_colors) <- anchor_ids
    sample_colors <- c(sample_colors, anchor_colors)
  }

  # Type colors
  type_colors <- c("μ-vector" = "#E31A1C", "Anchor" = "#1F78B4")

  annotation_colors <- list(
    Sample = sample_colors,
    Type = type_colors
  )

  # MODIFIED PATH: Save to final_report directory
  plots_dir <- here::here("final_report")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  timepoint_clean <- gsub("[^A-Za-z0-9]", "", timepoint_name)
  filename <- paste0("fig7_unified_heatmap_", timepoint_clean, "_", run_name, ".png")
  filepath <- file.path(plots_dir, filename)

  # Count samples for subtitle from the actual data used for the plot
  n_samples <- length(unique(mu_data$sample_info$sample_id[!mu_data$sample_info$is_anchor & mu_data$sample_info$sample_id != "Combined"]))

  # Create heatmap using pheatmap
  tryCatch({
    png(filepath, width = 14, height = 12, units = "in", res = 300)

    pheatmap(
      dist_matrix,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      annotation_row = annotation_df,
      annotation_col = annotation_df,
      annotation_colors = annotation_colors,
      main = paste0("Component and Anchor Euclidean Distance Heatmap\n",
                    timepoint_name, " (", n_samples, " samples + anchors)"),
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 8,
      show_rownames = TRUE,
      show_colnames = TRUE,
      angle_col = 45,
      cellwidth = 12,
      cellheight = 12,
      color = colorRampPalette(c("blue", "white", "red"))(100)
    )

    dev.off()

    if (file.exists(filepath)) {
      cat("✅ Unified heatmap saved to:", filepath, "\n")
    } else {
      cat("❌ Failed to create heatmap file\n")
    }

  }, error = function(e) {
    if(dev.cur() != 1) dev.off()  # Ensure device is closed
    cat("❌ Error creating heatmap:", e$message, "\n")
    return(NULL)
  })

  return(filepath)
}

#' Generate unified heatmaps for all timepoints
#' @param dataset_name Name of the dataset
#' @param run_name Run identifier
#' @param use_pcs Number of PC dimensions to use
#' @return List of file paths for all generated heatmaps
generate_all_unified_condition_heatmaps <- function(dataset_name = "HIV_PBMC",
                                                   run_name = "6-20k",
                                                   use_pcs = 10) {

  cat("🔥 UNIFIED CONDITION HEATMAP GENERATION\n")
  cat("======================================\n\n")

  # Load sample metadata
  sample_metadata <- load_sample_metadata_unified(dataset_name, run_name)

  # Get unique timepoints
  timepoints <- unique(sample_metadata$timepoint_clean)
  timepoints <- timepoints[!is.na(timepoints)]

  cat("Generating unified heatmaps for", length(timepoints), "timepoints:\n")
  cat(paste(timepoints, collapse = ", "), "\n")
  cat("Using Euclidean distance with first", use_pcs, "PCs\n\n")

  # Generate heatmaps for each timepoint
  all_filepaths <- list()

  for (timepoint in timepoints) {
    cat("Processing timepoint:", timepoint, "\n")

    filepath <- create_unified_timepoint_heatmap(
      dataset_name = dataset_name,
      timepoint_name = timepoint,
      sample_metadata = sample_metadata,
      run_name = run_name,
      use_pcs = use_pcs
    )

    if (!is.null(filepath)) {
      all_filepaths[[timepoint]] <- filepath
    }

    cat("\n")
  }

  cat("Generated", length(all_filepaths), "unified condition heatmaps\n\n")

  # Create summary statistics
  summary_stats <- sample_metadata %>%
    group_by(timepoint_clean) %>%
    summarise(
      n_samples = n(),
      patients = paste(sort(unique(patient)), collapse = ", "),
      k_range = paste(min(k_components, na.rm = TRUE), "to", max(k_components, na.rm = TRUE)),
      avg_k = round(mean(k_components, na.rm = TRUE), 1),
      models = paste(unique(model_type), collapse = ", "),
      .groups = "drop"
    )

  cat("=== UNIFIED HEATMAP GENERATION SUMMARY ===\n")
  print(summary_stats)

  return(list(
    filepaths = all_filepaths,
    metadata = sample_metadata,
    summary = summary_stats,
    settings = list(
      dataset_name = dataset_name,
      run_name = run_name,
      use_pcs = use_pcs
    )
  ))
}
