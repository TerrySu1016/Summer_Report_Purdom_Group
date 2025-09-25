# generate_6months_highlighted_plot.r - ARCHIVED VERSION
# Generate 6 Months PCA plot with specific components highlighted in red.
# This version is modified for the self-contained summer report archive.

# Load necessary libraries
library(here)
library(dplyr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(Seurat)

# MODIFIED PATH: Source the utility script from the local `./code/utils` directory
source(here::here("code", "utils", "util_multi_patient_visualization.r"))

#' Create modified 6 Months plot with highlighted components
#' @param dataset_name Name of the dataset
#' @param run_name Run identifier
#' @param highlighted_components Vector of component identifiers to highlight in red
#' @return ggplot object
create_highlighted_6months_plot <- function(dataset_name, run_name = "6-20k",
                                           highlighted_components = NULL) {

  # Load sample metadata
  sample_metadata <- load_sample_metadata(dataset_name, run_name)

  # Filter for 6 Months timepoint
  timepoint_samples <- sample_metadata %>%
    filter(timepoint_clean == "6  Months")

  if (nrow(timepoint_samples) == 0) {
    stop("No samples found for 6 Months timepoint")
  }

  cat("Creating plot for 6 Months with", nrow(timepoint_samples), "samples\n")

  # Load data for all samples in 6 Months timepoint
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

      # Create component identifier for highlighting
      mu_data$component_id <- paste0(sample_id, "_", mu_data$component)

      # Add highlighting flag
      mu_data$is_highlighted <- mu_data$component_id %in% highlighted_components

      all_mu_data <- rbind(all_mu_data, mu_data)
    }
  }

  if (nrow(all_pc_data) == 0) {
    stop("No PC data loaded for 6 Months timepoint")
  }

  # Get cell type colors
  cell_type_colors <- get_cell_type_colors(dataset_name)

  # Split mu data into highlighted and normal
  mu_normal <- all_mu_data[!all_mu_data$is_highlighted, ]
  mu_highlighted <- all_mu_data[all_mu_data$is_highlighted, ]

  cat("  Highlighting", nrow(mu_highlighted), "components in red\n")
  cat("  Normal components:", nrow(mu_normal), "\n")

  # Create the plot
  p <- ggplot() +
    # Cell data as background points
    geom_point(data = all_pc_data,
               aes(x = PC_1, y = PC_2, color = cell_type),
               alpha = 0.6, size = 0.5) +
    scale_color_manual(values = cell_type_colors, name = "Cell Type") +
    # Normal GMM mu vectors as black stars
    geom_point(data = mu_normal,
               aes(x = PC_1, y = PC_2),
               color = "black", shape = 8, size = 3, stroke = 1.5) +
    # Highlighted GMM mu vectors as red stars
    geom_point(data = mu_highlighted,
               aes(x = PC_1, y = PC_2),
               color = "red", shape = 8, size = 4, stroke = 2) +
    # Facet by patient
    facet_wrap(~ paste0(patient, " (K=", k_components, ")"),
               scales = "free", ncol = 2) +
    labs(
      title = "Multi-Patient PCA Visualization: 6 Months",
      subtitle = "Cells colored by type, black * = GMM components, red * = highlighted components",
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

# Set parameters
dataset_name <- "HIV_PBMC"
run_name <- "6-20k"

# Define components to highlight in red
# These components correspond to the "ideal" B-cell block identified in the analysis.
highlighted_components <- c(
  "P1_6 Months_C4",
  "P2_6 Months_C1",
  "P4_6 Months_C1"
)

cat("=== GENERATING 6 MONTHS PCA PLOT WITH HIGHLIGHTED COMPONENTS ===\n")
cat("Dataset:", dataset_name, "\n")
cat("Run:", run_name, "\n")
cat("Highlighted components:", paste(highlighted_components, collapse = ", "), "\n\n")

# Generate the plot
plot_obj <- create_highlighted_6months_plot(
  dataset_name = dataset_name,
  run_name = run_name,
  highlighted_components = highlighted_components
)

# MODIFIED PATH: Save the output to the `final_report` directory
output_dir <- here::here("final_report")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# MODIFIED FILENAME: Rename to match figure number in the report
plot_file <- file.path(output_dir, paste0("fig6_highlighted_6months_pca_", run_name, ".png"))
ggsave(plot_file, plot_obj, width = 12, height = 8, dpi = 300)

cat("Plot saved to:", plot_file, "\n")

# Display the plot
print(plot_obj)

cat("\nDone!\n")
