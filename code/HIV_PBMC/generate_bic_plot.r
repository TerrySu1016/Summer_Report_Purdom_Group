# Generate BIC Plot - ARCHIVED VERSION
# Creates BIC (Bayesian Information Criterion) plots for model selection.

here::i_am("code/HIV_PBMC/generate_bic_plot.r")

source(here::here("code", "utils", "util_bic_visualization.r"))

# Generate BIC analysis plot
bic_results <- quick_bic_analysis(
  dataset_name = "HIV_PBMC",
  run_name = "6-20k"
)
