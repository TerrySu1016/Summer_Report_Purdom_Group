# Unified Condition Heatmaps (Enhanced Style) - ARCHIVED VERSION
# Creates a single heatmap per timepoint containing all patients.
# This version is modified for the self-contained summer report archive.

here::i_am("code/HIV_PBMC/unified_condition_heatmaps.r")

# MODIFIED PATH: Source the utility script from the local `./code/utils` directory
source(here::here("code", "utils", "util_unified_condition_heatmaps.r"))

cat("🔥 STEP 15: UNIFIED CONDITION HEATMAPS (ENHANCED STYLE)\n")
cat("======================================================\n\n")

# This script is being run for a specific timepoint as part of the report.
# We will generate the heatmap only for the "6 Months" timepoint to match the report's narrative.
TARGET_TIMEPOINT <- "6 Months"

cat("🎯 CONFIGURATION:\n")
cat("• Distance metric: Euclidean (first 10 PCs)\n")
cat("• Timepoint focus:", TARGET_TIMEPOINT, "\n")
cat("• Output: High-resolution PNG file to `final_report` directory\n\n")

cat("🚀 RUNNING UNIFIED HEATMAP GENERATION FOR:", TARGET_TIMEPOINT, "\n")
cat("=============================================================\n")

# Load necessary metadata
sample_metadata <- load_sample_metadata_unified(dataset_name = "HIV_PBMC", run_name = "6-20k")

# Generate a single timepoint unified heatmap
# The output path is now controlled by the utility script to be `final_report/`
heatmap_filepath <- create_unified_timepoint_heatmap(
  dataset_name = "HIV_PBMC",
  timepoint_name = TARGET_TIMEPOINT,
  sample_metadata = sample_metadata,
  run_name = "6-20k",
  use_pcs = 10
)

cat("\n✅ UNIFIED HEATMAP GENERATION COMPLETED\n")
cat("=========================================\n\n")

if (!is.null(heatmap_filepath) && file.exists(heatmap_filepath)) {
  cat("📁 Generated Heatmap File:\n")
  cat("   •", TARGET_TIMEPOINT, ":", basename(heatmap_filepath), "\n")
} else {
  cat("❌ No unified heatmap was generated for", TARGET_TIMEPOINT, "\n")
  cat("   Check data availability and try again.\n")
}

cat("\n=== STEP 15 COMPLETION ===\n")
