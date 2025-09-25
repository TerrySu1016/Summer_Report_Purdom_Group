# util_visualize_gloscope.r - ARCHIVED VERSION
# Utility for creating MDS and heatmap plots from GloScope distance matrices.
# This version is modified for the self-contained summer report archive.
# Key change: All file paths now point to a local `./data` directory.

here::i_am("code/utils/util_visualize_gloscope.r")

# Load required libraries
library(GloScope)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)

#' Visualize GloScope Divergence Results
#' @param dataset_name Name of the dataset
#' @param methods Vector of method names (used in id_str when running gloscope)
#' @param run_name Name of the run (used in file naming)
#' @param metadata_columns Vector of column names to use from metadata (first is sample_id)
#' @param patient_id_col Column name for patient ID, for color coding
#' @param time_point_col Column name for time point, for shape coding
visualize_gloscope_results <- function(dataset_name, methods, run_name, metadata_columns,
                                     patient_id_col = NULL, time_point_col = NULL) {
  # Load Seurat object
  name_str <- "default_seurat"
  # MODIFIED PATH: Point to the local data directory
  seurat_object <- readRDS(here::here("data", "processed_data", paste0(name_str, ".rds")))
  message("Seurat object loaded successfully.")

  # Create sample metadata dataframe
  sample_id_col <- metadata_columns[1]
  all_req_cols <- unique(c(metadata_columns, patient_id_col, time_point_col))

  # Deduplicate metadata to have one row per sample
  all_metadata <- seurat_object@meta.data
  unique_metadata <- all_metadata[!duplicated(all_metadata[[sample_id_col]]), ]
  metadata_df <- unique_metadata[, all_req_cols]
  rownames(metadata_df) <- metadata_df[[sample_id_col]]

  message("Sample metadata dataframe created successfully.")

  # MODIFIED PATH: Save to final_report directory
  plots_dir <- here::here("final_report")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # Process each method's distance matrix
  for (method in methods) {
    # Load the saved distance matrix
    # MODIFIED PATH: Point to the local data directory
    dist_mat_file <- here::here("data", "gloscope_results",
                                  paste0(method, "_", run_name, "_distMat.rds"))
    if(!file.exists(dist_mat_file)) {
        stop("Distance matrix not found: ", dist_mat_file)
    }
    dist_mat <- readRDS(dist_mat_file)

    # 1. Generate MDS plot
    mds_result <- plotMDS(dist_mat = dist_mat,
                         metadata_df = metadata_df,
                         sample_id = sample_id_col,
                         color_by = patient_id_col,
                         shape_by = time_point_col,
                         k = 2)

    # Customize MDS plot
    mds_plot <- mds_result$plot +
      ggtitle(paste0("MDS plot of GloScope divergences (", method, ")")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    # Add manual color scale for patients
    if (!is.null(patient_id_col)) {
        color_values <- c("P1" = "#1f77b4", "P2" = "#ff7f0e", "P3" = "#2ca02c", "P4" = "#d62728")
        mds_plot <- mds_plot +
            scale_color_manual(name = "Patient", values = color_values)
    }

    # Add manual shape scale for time points
    if (!is.null(time_point_col)) {
        shape_values <- c(
            "0 Weeks" = 16,        # filled circle
            "1 Week" = 17,         # filled triangle
            "1 Year" = 15,         # filled square
            "2 Weeks" = 3,         # plus
            "3 Weeks" = 4,         # cross
            "4 Weeks" = 8,         # asterisk
            "6 Months" = 18,       # filled diamond
            "Pre-Infection" = 5    # open diamond
        )
        mds_plot <- mds_plot +
            scale_shape_manual(name = "Time Point", values = shape_values)
    }

    # MODIFIED FILENAME: Save as PNG with figure number
    ggsave(filename = file.path(plots_dir, paste0("fig4_mds_", method, "_", run_name, ".png")),
           plot = mds_plot,
           width = 7,
           height = 6,
           dpi = 300)
           
    message("MDS plot saved to ", plots_dir)

    # 2. Generate heatmap
    group_id_col <- patient_id_col %||% metadata_columns[2]
    ordered_samples <- metadata_df[order(metadata_df[[group_id_col]]), sample_id_col]
    dist_mat_ordered <- dist_mat[ordered_samples, ordered_samples]

    # Generate heatmap without patient annotations
    pheatmap(dist_mat_ordered,
            color = viridis(100),
            main = paste0("GloScope divergences (", method, ")"),
            show_rownames = TRUE,
            show_colnames = TRUE,
            filename = file.path(plots_dir, paste0("heatmap_", method, "_", run_name, ".png")),
            width = 8,
            height = 7)
    message("Heatmap saved to ", plots_dir)
  }
  cat("Visualization complete. Results saved to:", plots_dir, "\n")
}
