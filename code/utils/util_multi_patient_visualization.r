# util_multi_patient_visualization.r - ARCHIVED VERSION
# Creates timepoint-grouped faceted plots showing all patients per condition.
# This version is modified for the self-contained summer report archive.
# Key change: All file paths now point to a local `./data` directory.

here::i_am("code/utils/util_multi_patient_visualization.r")

library(here)
library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
library(RColorBrewer)
library(Seurat)

#' Get global cell types and create color mapping
#' @param dataset_name Name of the dataset
#' @return Named vector of colors for each cell type
get_cell_type_colors <- function(dataset_name = "HIV_PBMC") {

  # MODIFIED PATH: Point to the local data directory
  seurat_file <- here::here("data", "processed_data", "default_seurat.rds")

  if (!file.exists(seurat_file)) {
    stop("Seurat object not found: ", seurat_file)
  }

  seurat_obj <- readRDS(seurat_file)
  all_cell_types <- sort(unique(seurat_obj@meta.data$cell_type))

  # Create consistent color mapping
  predefined_colors <- c(
    "B cell" = "#FFD700",              # Bright gold
    "T cell" = "#8A2BE2",              # Blue violet
    "Monocyte" = "#FF69B4",            # Hot pink
    "NK cell" = "#00CED1",             # Dark turquoise
    "Dendritic cell" = "#FF4500",      # Orange red
    "Plasmablast" = "#32CD32",         # Lime green
    "Platelet" = "#DC143C"             # Crimson
  )

  # Assign colors to all cell types
  cell_type_colors <- rep(NA, length(all_cell_types))
  names(cell_type_colors) <- all_cell_types

  # Use predefined colors where available
  for (cell_type in names(predefined_colors)) {
    if (cell_type %in% all_cell_types) {
      cell_type_colors[cell_type] <- predefined_colors[cell_type]
    }
  }

  # Generate additional colors for any remaining cell types
  remaining_types <- all_cell_types[is.na(cell_type_colors)]
  if (length(remaining_types) > 0) {
    additional_colors <- rainbow(length(remaining_types))
    cell_type_colors[remaining_types] <- additional_colors
  }

  return(cell_type_colors)
}

#' Load available samples and group by timepoint
#' @param dataset_name Name of the dataset
#' @param run_name Run identifier
#' @return Data frame with sample metadata grouped by timepoint
load_sample_metadata <- function(dataset_name = "HIV_PBMC", run_name = "6-20k") {

  # MODIFIED PATH: Point to the local data directory with run-specific subdirectory
  data_dir <- here::here("data", "gloscope_results", run_name)

  # Look for GMM parameter files to identify available samples
  param_files <- list.files(data_dir, pattern = paste0("gmm_params_.*_defaultGMM_", run_name, "\\.rds"), full.names = TRUE)

  if (length(param_files) == 0) {
    stop("No GMM parameter files found in ", data_dir)
  }

  # Extract sample IDs from file names
  sample_ids <- gsub(paste0(".*gmm_params_(.*)_defaultGMM_", run_name, "\\.rds"), "\\1", basename(param_files))

  # Parse sample IDs to extract patient and timepoint information
  sample_info <- data.frame(
    sample_id = sample_ids,
    stringsAsFactors = FALSE
  )

  # Extract patient and timepoint from sample ID
  sample_info$patient <- gsub("_.*", "", sample_info$sample_id)
  sample_info$timepoint <- gsub(".*_", "", sample_info$sample_id)

  # Standardize timepoint names
  sample_info$timepoint_clean <- sample_info$timepoint
  sample_info$timepoint_clean <- gsub("Weeks?", " Weeks", sample_info$timepoint_clean)
  sample_info$timepoint_clean <- gsub("Months?", " Months", sample_info$timepoint_clean)
  sample_info$timepoint_clean <- gsub("Years?", " Year", sample_info$timepoint_clean)
  sample_info$timepoint_clean <- gsub("PreInfection", "Pre-Infection", sample_info$timepoint_clean)

  # Check if parameter files exist and load K values
  sample_info$k_components <- NA
  sample_info$model_type <- NA

  for (i in 1:nrow(sample_info)) {
    # MODIFICATION: Use data_dir
    param_file <- file.path(data_dir, paste0("gmm_params_", sample_info$sample_id[i], "_defaultGMM_", run_name, ".rds"))

    if (file.exists(param_file)) {
      tryCatch({
        params <- readRDS(param_file)
        sample_info$k_components[i] <- ncol(params$mu)  # Number of components
        sample_info$model_type[i] <- params$modelName %||% "Unknown"
      }, error = function(e) {
        cat("Warning: Could not load parameters for", sample_info$sample_id[i], "\n")
      })
    }
  }

  # Sort by timepoint and patient
  sample_info <- sample_info %>%
    arrange(timepoint_clean, patient)

  cat("Loaded metadata for", nrow(sample_info), "samples\n")
  cat("Timepoints:", paste(unique(sample_info$timepoint_clean), collapse = ", "), "\n")
  cat("K values range:", min(sample_info$k_components, na.rm = TRUE), "to", max(sample_info$k_components, na.rm = TRUE), "\n\n")

  return(sample_info)
}

#' Load PC data for a sample from Seurat object
#' @param dataset_name Name of the dataset
#' @param sample_id Sample identifier
#' @param run_name Run identifier
#' @return Data frame with PC coordinates and cell type information
load_sample_pca_data <- function(dataset_name, sample_id, run_name = "6-20k") {

  # MODIFIED PATH: Point to the local data directory with run-specific subdirectory
  pc_file <- here::here("data", "gloscope_results", run_name,
                        paste0("pc_coords_", sample_id, "_", run_name, ".csv"))

  if (file.exists(pc_file)) {
    pc_data <- read_csv(pc_file, show_col_types = FALSE)
    pc_data$sample_id <- sample_id
    return(pc_data)
  }

  # If CSV doesn't exist, extract from Seurat object
  # MODIFIED PATH: Point to the local data directory
  seurat_file <- here::here("data", "processed_data", "default_seurat.rds")

  if (!file.exists(seurat_file)) {
    warning("Neither PC file nor Seurat object found for sample ", sample_id)
    return(NULL)
  }

  tryCatch({
    # Load Seurat object
    seurat_obj <- readRDS(seurat_file)

    # Filter cells for this sample
    sample_cells <- seurat_obj@meta.data$sample == sample_id

    if (sum(sample_cells) == 0) {
      warning("No cells found for sample ", sample_id, " in Seurat object")
      return(NULL)
    }

    # Extract PC coordinates
    pca_embeddings <- seurat_obj@reductions$pca@cell.embeddings[sample_cells, 1:10]

    # Extract cell type information
    cell_types <- seurat_obj@meta.data$cell_type[sample_cells]

    # Create data frame
    pc_data <- as.data.frame(pca_embeddings)
    pc_data$cell_type <- cell_types
    pc_data$sample_id <- sample_id

    # Rename PC columns to match expected format
    colnames(pc_data)[1:10] <- paste0("PC_", 1:10)

    cat("  Extracted", nrow(pc_data), "cells for", sample_id, "from Seurat object\n")

    return(pc_data)

  }, error = function(e) {
    warning("Error extracting PC data for ", sample_id, ": ", e$message)
    return(NULL)
  })
}

#' Load GMM mu vectors for a sample
#' @param dataset_name Name of the dataset
#' @param sample_id Sample identifier
#' @param run_name Run identifier
#' @return Data frame with mu vector coordinates
load_sample_mu_vectors <- function(dataset_name, sample_id, run_name = "6-20k") {

  # MODIFIED PATH: Point to the local data directory with run-specific subdirectory
  param_file <- here::here("data", "gloscope_results", run_name,
                           paste0("gmm_params_", sample_id, "_defaultGMM_", run_name, ".rds"))

  if (!file.exists(param_file)) {
    warning("Parameter file not found for sample ", sample_id)
    return(NULL)
  }

      tryCatch({
      params <- readRDS(param_file)

      if (!is.null(params$mu)) {
        mu_matrix <- params$mu

        # Convert to data frame
        if (is.matrix(mu_matrix)) {
          mu_df <- as.data.frame(t(mu_matrix))  # Transpose so each row is a component

          # Set PC column names first (before adding other columns)
          n_pcs <- nrow(mu_matrix)  # Number of PCs (rows in original matrix)
          pc_cols <- paste0("PC_", 1:n_pcs)
          colnames(mu_df) <- pc_cols

          # Now add metadata columns
          mu_df$component <- paste0("C", 1:nrow(mu_df))
          mu_df$sample_id <- sample_id

          return(mu_df)
        }
      }

      return(NULL)

    }, error = function(e) {
      warning("Error loading mu vectors for ", sample_id, ": ", e$message)
      return(NULL)
    })
}

#' Create faceted PCA plot for all patients in a timepoint
#' @param dataset_name Name of the dataset
#' @param timepoint_name Name of the timepoint (e.g., "0 Weeks", "1 Week")
#' @param sample_metadata Data frame with sample metadata
#' @param run_name Run identifier
#' @param cell_type_colors Named vector of colors for cell types
#' @return ggplot object
create_timepoint_faceted_plot <- function(dataset_name, timepoint_name, sample_metadata, run_name = "6-20k", cell_type_colors = NULL) {

  # Get samples for this timepoint
  timepoint_samples <- sample_metadata %>%
    filter(timepoint_clean == timepoint_name)

  if (nrow(timepoint_samples) == 0) {
    warning("No samples found for timepoint: ", timepoint_name)
    return(NULL)
  }

  cat("Creating plot for", timepoint_name, "with", nrow(timepoint_samples), "samples\n")

  # Load data for all samples in this timepoint
  all_pc_data <- data.frame()
  all_mu_data <- data.frame()

  for (i in 1:nrow(timepoint_samples)) {
    sample_id <- timepoint_samples$sample_id[i]
    patient <- timepoint_samples$patient[i]
    k_comp <- timepoint_samples$k_components[i]

    cat("  Loading", sample_id, "(", patient, ", K =", k_comp, ")\n")

    # Load PC data
    pc_data <- load_sample_pca_data(dataset_name, sample_id, run_name)
    if (!is.null(pc_data)) {
      pc_data$patient <- patient
      pc_data$k_components <- k_comp
      all_pc_data <- rbind(all_pc_data, pc_data)
    }

    # Load mu vectors
    mu_data <- load_sample_mu_vectors(dataset_name, sample_id, run_name)
    if (!is.null(mu_data)) {
      mu_data$patient <- patient
      mu_data$k_components <- k_comp
      all_mu_data <- rbind(all_mu_data, mu_data)
    }
  }

  if (nrow(all_pc_data) == 0) {
    warning("No PC data loaded for timepoint: ", timepoint_name)
    return(NULL)
  }

  # Get cell type colors if not provided
  if (is.null(cell_type_colors)) {
    cell_type_colors <- get_cell_type_colors(dataset_name)
  }

  # Create the plot
  p <- ggplot() +
    # Cell data as background points
    geom_point(data = all_pc_data,
               aes(x = PC_1, y = PC_2, color = cell_type),
               alpha = 0.6, size = 0.5) +
    scale_color_manual(values = cell_type_colors, name = "Cell Type") +
    # GMM mu vectors as black stars (no color mapping needed)
    geom_point(data = all_mu_data,
               aes(x = PC_1, y = PC_2),
               color = "black", shape = 8, size = 3, stroke = 1.5) +
    # Facet by patient
    facet_wrap(~ paste0(patient, " (K=", k_components, ")"),
               scales = "free", ncol = 2) +
    labs(
      title = paste("Multi-Patient PCA Visualization:", timepoint_name),
      subtitle = "Cells colored by type, black * = GMM components (mu vectors)",
      x = "PC_1",
      y = "PC_2"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10)
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

  return(p)
}

#' Generate all timepoint plots
#' @param dataset_name Name of the dataset
#' @param run_name Run identifier
#' @param save_plots Whether to save plots to files
#' @param output_dir Directory to save plots (if NULL, uses default)
#' @return List of ggplot objects
generate_all_timepoint_plots <- function(dataset_name = "HIV_PBMC",
                                        run_name = "6-20k",
                                        save_plots = TRUE,
                                        output_dir = NULL) {

  cat("🎯 MULTI-PATIENT TIMEPOINT VISUALIZATION\n")
  cat("=======================================\n\n")

  # Load sample metadata
  sample_metadata <- load_sample_metadata(dataset_name, run_name)

  # Get cell type colors
  cell_type_colors <- get_cell_type_colors(dataset_name)

  # Get unique timepoints
  timepoints <- unique(sample_metadata$timepoint_clean)
  timepoints <- timepoints[!is.na(timepoints)]

  cat("Generating plots for", length(timepoints), "timepoints:\n")
  cat(paste(timepoints, collapse = ", "), "\n\n")

  # Set up output directory
  if (is.null(output_dir)) {
    # MODIFIED PATH: Point to the final_report directory in the archive
    output_dir <- here::here("final_report")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate plots for each timepoint
  all_plots <- list()

  for (timepoint in timepoints) {
    cat("Processing timepoint:", timepoint, "\n")

    plot_obj <- create_timepoint_faceted_plot(
      dataset_name = dataset_name,
      timepoint_name = timepoint,
      sample_metadata = sample_metadata,
      run_name = run_name,
      cell_type_colors = cell_type_colors
    )

    if (!is.null(plot_obj)) {
      all_plots[[timepoint]] <- plot_obj

      # Save plot if requested
      if (save_plots) {
        # Clean timepoint name for filename
        clean_timepoint <- gsub("[^A-Za-z0-9]", "_", timepoint)
        plot_file <- file.path(output_dir, paste0("multi_patient_pca_", clean_timepoint, "_", run_name, ".png"))

        ggsave(plot_file, plot_obj, width = 12, height = 8, dpi = 300)
        cat("  Plot saved to:", plot_file, "\n")
      }

      # Display plot
      print(plot_obj)
    }

    cat("\n")
  }

  cat("Generated", length(all_plots), "timepoint plots\n\n")

  # Create summary table
  summary_table <- sample_metadata %>%
    group_by(timepoint_clean) %>%
    summarise(
      n_samples = n(),
      patients = paste(sort(unique(patient)), collapse = ", "),
      k_range = paste(min(k_components, na.rm = TRUE), "to", max(k_components, na.rm = TRUE)),
      models = paste(unique(model_type), collapse = ", "),
      .groups = "drop"
    )

  cat("=== TIMEPOINT SUMMARY ===\n")
  print(summary_table)

  return(list(
    plots = all_plots,
    metadata = sample_metadata,
    summary = summary_table
  ))
}
