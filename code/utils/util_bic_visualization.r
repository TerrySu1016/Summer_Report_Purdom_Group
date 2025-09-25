# util_bic_visualization.r - ARCHIVED VERSION
# An enhanced BIC (Bayesian Information Criterion) visualization utility.
# This version is modified for the self-contained summer report archive.
# Key changes:
#   - All file paths now point to a local `./data` directory.
#   - All comments and console messages are in English.

here::i_am("code/utils/util_bic_visualization.r")

suppressMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

#' Extract BIC data from fitted models (Working Version)
#' @param dataset_name Name of the dataset
#' @param run_name Run identifier
#' @param k_range Range of K values to include
#' @param models Vector of model names to include
#' @param max_samples Maximum number of samples to process
#' @return Data frame with BIC values
extract_bic <- function(dataset_name = "HIV_PBMC",
                               run_name = "6-20k",
                               k_range = 6:20,
                               models = c("EVE", "VVE", "VVV"),
                               max_samples = 10) {

  cat("📊 Extracting BIC data from fitted models\n")
  cat("=========================================\n")

  # MODIFIED PATH: Point to the local data directory
  fitted_models_file <- here::here("data", "gloscope_results",
                                   paste0("gmm_fitted_models_", run_name, ".rds"))

  if (!file.exists(fitted_models_file)) {
    stop("❌ Fitted models file not found: ", fitted_models_file)
  }

  fitted_models <- readRDS(fitted_models_file)
  cat("✅ Loaded", length(fitted_models), "fitted models\n")

  # Limit number of samples for performance and select better samples
  model_keys <- names(fitted_models)
  
  # Prefer P4_6 Months over P3_6 Months for better data quality
  if ("P3_6 Months_defaultGMM" %in% model_keys && "P4_6 Months_defaultGMM" %in% model_keys) {
    # Remove P3_6 Months and keep P4_6 Months
    model_keys <- model_keys[model_keys != "P3_6 Months_defaultGMM"]
    cat("📊 Using P4_6 Months instead of P3_6 Months (better data quality)\n")
  }
  
  if (length(model_keys) > max_samples) {
    model_keys <- model_keys[1:max_samples]
    cat("📌 Limiting processing to", max_samples, "samples for performance\n")
  }

  # Extract BIC data
  bic_data <- data.frame()

  for (model_key in model_keys) {
    model <- fitted_models[[model_key]]

    if (!is.null(model$BIC)) {
      bic_matrix <- model$BIC

      # Extract sample info from model key
      sample_parts <- strsplit(model_key, "_")[[1]]
      sample_id <- paste(sample_parts[1:(length(sample_parts)-1)], collapse="_")

      cat("  Processing sample:", sample_id, "\n")

      # Get available models
      available_models <- intersect(models, colnames(bic_matrix))

      # Extract data for specified K range
      for (k in k_range) {
        if (as.character(k) %in% rownames(bic_matrix)) {
          for (model_name in available_models) {
            bic_value <- bic_matrix[as.character(k), model_name]

            # More robust checking for valid BIC values
            if (!is.na(bic_value) && is.finite(bic_value) && !is.null(bic_value)) {
              bic_data <- rbind(bic_data, data.frame(
                sample_id = sample_id,
                model = model_name,
                k = k,
                bic = bic_value,
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
  }

  cat("✅ Extraction complete! Total of", nrow(bic_data), "BIC values\n")
  cat("   Samples:", length(unique(bic_data$sample_id)), "\n")
  cat("   Models:", paste(unique(bic_data$model), collapse=", "), "\n")
  cat("   K range:", min(bic_data$k), "-", max(bic_data$k), "\n")
  
  # Add diagnostic information about data coverage
  if (nrow(bic_data) > 0) {
    data_coverage <- bic_data %>%
      group_by(sample_id, k) %>%
      summarise(n_models = n(), .groups = "drop") %>%
      group_by(sample_id) %>%
      summarise(
        k_range = paste(min(k), "-", max(k)),
        avg_models_per_k = round(mean(n_models), 1),
        .groups = "drop"
      )
    
    cat("\n📊 Data coverage by sample:\n")
    for (i in 1:nrow(data_coverage)) {
      cat(sprintf("   %s: K=%s, avg %.1f models per K\n", 
                  data_coverage$sample_id[i], 
                  data_coverage$k_range[i], 
                  data_coverage$avg_models_per_k[i]))
    }
  }
  cat("\n")

  return(bic_data)
}

#' Create enhanced BIC plot (Working Version)
#' @param bic_data Data frame from extract_bic()
#' @param highlight_best Whether to highlight optimal K points
#' @param color_palette Color palette to use
#' @param facet_by_sample Whether to create separate panels for each sample
#' @return ggplot object
create_bic_plot <- function(bic_data,
                                   highlight_best = TRUE,
                                   color_palette = "viridis",
                                   facet_by_sample = TRUE) {

  cat("📊 Creating enhanced BIC plot\n")
  cat("============================\n")

  # Determine the full K range from the data for setting limits
  k_limits <- range(bic_data$k, na.rm = TRUE)

  # Prepare data with optimal K annotations - handle missing data gracefully
  bic_enhanced <- bic_data %>%
    # Remove any remaining invalid values
    filter(!is.na(bic) & is.finite(bic)) %>%
    group_by(sample_id, model) %>%
    mutate(
      is_best_k = if(n() > 0) bic == max(bic, na.rm = TRUE) else FALSE,
      best_k_this_model = if(n() > 0) k[which.max(bic)] else NA
    ) %>%
    ungroup() %>%
    group_by(sample_id) %>%
    mutate(
      global_best = if(n() > 0) bic == max(bic, na.rm = TRUE) else FALSE,
      best_model_overall = if(n() > 0) model[which.max(bic)] else NA,
      best_k_overall = if(n() > 0) k[which.max(bic)] else NA
    ) %>%
    ungroup()

  # Create base plot with better handling of missing data
  p <- ggplot(bic_enhanced, aes(x = k, y = bic, color = model)) +
    geom_line(linewidth = 1.3, alpha = 0.8, na.rm = TRUE) +
    geom_point(size = 2.5, alpha = 0.9, na.rm = TRUE) +
    # Force the x-axis to show the complete K range even if some data is missing
    scale_x_continuous(limits = c(6, 20), breaks = seq(6, 20, 2))

  # Apply color scheme
  if (color_palette == "viridis") {
    p <- p + scale_color_viridis_d(name = "Model")
  } else if (color_palette == "Set2") {
    p <- p + scale_color_brewer(type = "qual", palette = "Set2",
                               name = "Model")
  } else {
    # Default colors
    colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A")
    p <- p + scale_color_manual(values = colors,
                               name = "Model")
  }

  # Highlight optimal points
  if (highlight_best) {
    # Highlight best K for each model
    best_k_points <- bic_enhanced %>% filter(is_best_k)
    p <- p +
      geom_point(data = best_k_points,
                 aes(x = k, y = bic, color = model),
                 size = 4, shape = 21, fill = "white", stroke = 2) +
      geom_text(data = best_k_points,
                aes(x = k, y = bic, label = paste0("K=", k)),
                vjust = -1.2, hjust = 0.5, size = 3.5, fontface = "bold",
                color = "black")

    # Highlight global best points
    global_best_points <- bic_enhanced %>% filter(global_best)
    if (nrow(global_best_points) > 0) {
      p <- p +
        geom_point(data = global_best_points,
                   aes(x = k, y = bic),
                   size = 6, shape = 23, fill = "gold", color = "black", stroke = 2)
    }
  }

  # Add faceting if requested
  if (facet_by_sample && length(unique(bic_enhanced$sample_id)) > 1) {
    p <- p + facet_wrap(~ sample_id, scales = "free_y", ncol = 2)
  }

  # Add labels and theme
    p <- p + labs(
      title = "BIC Curves: Covariance Structure Model Comparison",
      subtitle = if(highlight_best) "◇ = Global optimal, ○ = Model optimal K (Higher BIC = Better)" else "Higher BIC = Better fit",
      x = "Number of Components (K)",
      y = "BIC Value",
      caption = paste("Analysis based on", length(unique(bic_enhanced$sample_id)), "samples")
    )

  p <- p + theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "darkgray"),
      plot.caption = element_text(size = 10, color = "darkgray", hjust = 1),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90")
    )

  return(p)
}


#' Generate summary of optimal models
#' @param bic_data Data frame with BIC values
#' @return Data frame with optimal model summary
summarize_optimal_models <- function(bic_data) {

  summary_stats <- bic_data %>%
    group_by(sample_id) %>%
    slice_max(bic, n = 1) %>%
    ungroup() %>%
    select(sample_id, model, k, bic) %>%
    arrange(desc(bic))

  cat("🔍 Optimal Model Summary\n")
  cat("=========================\n")
  print(summary_stats)
  cat("\n")

  # Model frequency
  model_freq <- summary_stats %>%
    count(model, sort = TRUE) %>%
    mutate(percentage = round(n/sum(n)*100, 1))

  cat("📊 Optimal Model Distribution:\n")
  for (i in 1:nrow(model_freq)) {
    cat(sprintf("  %s: %d times (%.1f%%)\n",
                model_freq$model[i],
                model_freq$n[i],
                model_freq$percentage[i]))
  }
  cat("\n")

  return(summary_stats)
}

#' Complete enhanced BIC analysis (Working Version)
#' @param dataset_name Name of the dataset
#' @param run_name Run identifier
#' @param k_range Range of K values to test
#' @param models Vector of model names to test
#' @param max_samples Maximum number of samples to process
#' @param save_results Whether to save results
#' @param color_palette Color palette
#' @return List with plot, data, and summary
run_bic_analysis <- function(dataset_name = "HIV_PBMC",
                                    run_name = "6-20k",
                                    k_range = 6:20,
                                    models = c("EVE", "VVE", "VVV"),
                                    max_samples = 8,
                                    save_results = TRUE,
                                    color_palette = "viridis") {

  cat("🎨 Enhanced BIC Analysis\n")
  cat("=========================================\n\n")

  # Extract BIC data
  bic_data <- extract_bic(
    dataset_name = dataset_name,
    run_name = run_name,
    k_range = k_range,
    models = models,
    max_samples = max_samples
  )

  if (nrow(bic_data) == 0) {
    stop("❌ No BIC data found.")
  }

  # Create plot
  enhanced_plot <- create_bic_plot(
    bic_data = bic_data,
    highlight_best = TRUE,
    color_palette = color_palette,
    facet_by_sample = TRUE
  )

  # Generate summary
  summary_stats <- summarize_optimal_models(bic_data)

  # Save results if requested
  if (save_results) {
    # MODIFIED PATH: Save to final_report directory
    output_dir <- here::here("final_report")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Save plot
    plot_file <- file.path(output_dir, paste0("fig5_enhanced_bic_curves_", run_name, ".png"))
    ggsave(plot_file, enhanced_plot, width = 14, height = 10, dpi = 300)
    cat("📁 Plot saved to:", plot_file, "\n")

    # Save data
    data_file <- file.path(output_dir, paste0("bic_analysis_data_", run_name, ".csv"))
    write_csv(bic_data, data_file)
    cat("📊 Data saved to:", data_file, "\n")

    # Save summary
    summary_file <- file.path(output_dir, paste0("bic_analysis_summary_", run_name, ".csv"))
    write_csv(summary_stats, summary_file)
    cat("📋 Summary saved to:", summary_file, "\n\n")
  }

  # Display plot
  print(enhanced_plot)

  return(list(
    plot = enhanced_plot,
    data = bic_data,
    summary = summary_stats
  ))
}

#' Intelligent BIC analysis with adaptive K range
#' @param dataset_name Name of the dataset
#' @param run_name Run identifier
#' @param adaptive_k_range Whether to automatically adjust K range based on data availability
#' @return Complete BIC analysis results
quick_bic_analysis <- function(dataset_name = "HIV_PBMC",
                              run_name = "6-20k",
                              adaptive_k_range = TRUE) {

  # If adaptive K range is enabled, first check data availability
  effective_k_range <- 6:20
  if (adaptive_k_range) {
    cat("🔍 Checking data availability for adaptive K range...\n")
    
    # Load fitted models to check available K values
    fitted_models_file <- here::here("data", "gloscope_results",
                                     paste0("gmm_fitted_models_", run_name, ".rds"))
    
    if (file.exists(fitted_models_file)) {
      fitted_models <- readRDS(fitted_models_file)
      
      # Check K value coverage across samples
      k_coverage <- data.frame()
      for (model_key in names(fitted_models)[1:min(6, length(fitted_models))]) {
        model <- fitted_models[[model_key]]
        if (!is.null(model$BIC)) {
          bic_matrix <- model$BIC
          available_k <- as.numeric(rownames(bic_matrix))
          
          # Count valid (non-NA) values for key models
          key_models <- c("EVE", "VVE", "VVV")
          for (k in available_k) {
            valid_count <- sum(!is.na(bic_matrix[as.character(k), 
                                              intersect(key_models, colnames(bic_matrix))]))
            k_coverage <- rbind(k_coverage, data.frame(
              sample = model_key,
              k = k,
              valid_models = valid_count
            ))
          }
        }
      }
      
      if (nrow(k_coverage) > 0) {
        # Find K range where most samples have at least 2 valid models
        k_summary <- k_coverage %>%
          group_by(k) %>%
          summarise(
            samples_with_2plus = sum(valid_models >= 2),
            avg_valid_models = mean(valid_models),
            .groups = "drop"
          )
        
        # Use K values where at least 80% of samples have 2+ valid models
        min_samples_threshold <- max(1, 0.8 * length(unique(k_coverage$sample)))
        suitable_k <- k_summary$k[k_summary$samples_with_2plus >= min_samples_threshold]
        
        if (length(suitable_k) > 0) {
          effective_k_range <- min(suitable_k):max(suitable_k)
          cat(sprintf("📊 Adaptive K range: %d-%d (based on data coverage)\n", 
                      min(effective_k_range), max(effective_k_range)))
        } else {
          cat("⚠️ Using default K range due to limited data coverage\n")
        }
      }
    }
  }

  results <- run_bic_analysis(
    dataset_name = dataset_name,
    run_name = run_name,
    k_range = effective_k_range,
    models = c("EVE", "VVE", "VVV"),
    max_samples = 6,  # Process 6 representative samples
    save_results = TRUE,
    color_palette = "viridis"
  )

  return(results)
}
