# generate_enhanced_mds.r - ARCHIVED VERSION
# Generates MDS plots from pre-computed GloScope distance matrices.
# This version is modified for the self-contained summer report archive.

# Load necessary libraries
library(here)

# MODIFIED PATH: Source the utility script from the local `./code/utils` directory
source(here::here("code", "utils", "util_visualize_gloscope.r"))

# Define common parameters
dataset_name <- "HIV_PBMC"
# 'defaultGMM' is the method name used when the distance matrix was saved.
methods <- c("defaultGMM")

# Define metadata columns to be used from the Seurat object.
# The first element must be the sample identifier.
metadata_cols <- c("sample", "patient", "TimePoint")

# Specify the columns for color and shape aesthetics in the plot
patient_col <- "patient"
timepoint_col <- "TimePoint"

# --- Generate plots for 6-20k (Expanded k-range) ---
# This corresponds to the final MDS plot used in the report (Figure 4, right panel)
message("--- Generating enhanced MDS plot for 6-20k (Expanded k-range) ---")
visualize_gloscope_results(
  dataset_name = dataset_name,
  methods = methods,
  run_name = "6-20k",
  metadata_columns = metadata_cols,
  patient_id_col = patient_col,
  time_point_col = timepoint_col
)
message("Plot for 6-20k generated successfully in 'final_report' folder.")

# The original report also showed a comparison to a run with a different k-range.
# Since this archive is self-contained with a single set of data, we only
# generate the plot for the primary run. If data for '1-9k' were present in the
# `./data` folder, the following code could be activated to generate it.

# --- (Optional) Generate plots for a second run if data is available ---
# message("--- Generating enhanced MDS plots for 1-9k (Default k-range) ---")
# visualize_gloscope_results(
#   dataset_name = dataset_name,
#   methods = methods,
#   run_name = "1-9k",
#   metadata_columns = metadata_cols,
#   patient_id_col = patient_col,
#   time_point_col = timepoint_col
# )
# message("Plot for 1-9k generated successfully.")

message("\nAll MDS plots have been generated.")
