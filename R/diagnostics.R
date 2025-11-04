#' Diagnose data and recommend appropriate model
#'
#' Analyzes count or PA data to provide model recommendations based on
#' data characteristics like sparsity, overdispersion, and zero-inflation.
#'
#' @param data_matrix Matrix or data.frame of count or PA data (features x samples)
#' @param data_type Character: "auto" (default), "count", or "PA". If "auto", will detect based on values.
#' @param metadata Optional data frame with sample metadata for additional context
#' @param make_plots Logical; if TRUE, generates diagnostic plots (default: TRUE)
#' @param n_features_plot Number of features to show in example plots (default: 12)
#'
#' @return A list with:
#' \describe{
#'   \item{data_type}{Detected or specified data type}
#'   \item{summary_stats}{Summary statistics for the dataset}
#'   \item{recommendations}{Model recommendations with rationale}
#'   \item{diagnostics}{Detailed diagnostic metrics}
#'   \item{plots}{List of diagnostic plots (if make_plots = TRUE)}
#' }
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point geom_abline geom_density theme_minimal labs facet_wrap
#' @importFrom dplyr %>% mutate
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' \dontrun{
#' # Diagnose count data
#' diagnosis <- diagnose_data(my_count_matrix)
#'
#' # View recommendations
#' diagnosis$recommendations
#'
#' # View plots
#' diagnosis$plots$zero_distribution
#' diagnosis$plots$mean_variance
#' }
diagnose_data <- function(data_matrix,
                         data_type = c("auto", "count", "PA"),
                         metadata = NULL,
                         make_plots = TRUE,
                         n_features_plot = 12) {

  data_type <- match.arg(data_type)

  # Convert to matrix if needed
  if (is.data.frame(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }

  # Auto-detect data type
  if (data_type == "auto") {
    max_val <- max(data_matrix, na.rm = TRUE)
    unique_vals <- length(unique(as.vector(data_matrix)))

    if (max_val <= 1 && unique_vals <= 3) {
      data_type <- "PA"
      message("Detected data type: Presence/Absence (binary)")
    } else {
      data_type <- "count"
      message("Detected data type: Count/Abundance")
    }
  }

  # Initialize results
  results <- list(
    data_type = data_type,
    summary_stats = list(),
    diagnostics = list(),
    recommendations = list(),
    plots = list()
  )

  # Calculate summary statistics
  n_features <- nrow(data_matrix)
  n_samples <- ncol(data_matrix)

  results$summary_stats <- list(
    n_features = n_features,
    n_samples = n_samples,
    total_observations = n_features * n_samples
  )

  # Data-type specific diagnostics
  if (data_type == "PA") {
    results <- diagnose_PA_data(data_matrix, results, make_plots, n_features_plot)
  } else {
    results <- diagnose_count_data(data_matrix, results, make_plots, n_features_plot)
  }

  # Add metadata summary if provided
  if (!is.null(metadata)) {
    results$summary_stats$n_metadata_vars <- ncol(metadata)
    results$summary_stats$metadata_samples <- nrow(metadata)
    results$summary_stats$sample_overlap <- sum(colnames(data_matrix) %in% rownames(metadata))
  }

  class(results) <- c("longitudinal_diagnosis", "list")
  return(results)
}

#' Diagnose PA data
#' @keywords internal
#' @noRd
diagnose_PA_data <- function(data_matrix, results, make_plots, n_features_plot) {

  # Calculate prevalence per feature
  prevalence <- rowMeans(data_matrix, na.rm = TRUE)

  results$diagnostics$prevalence <- list(
    mean = mean(prevalence),
    median = median(prevalence),
    min = min(prevalence),
    max = max(prevalence),
    very_rare = sum(prevalence < 0.05) / length(prevalence),
    rare = sum(prevalence < 0.10) / length(prevalence),
    common = sum(prevalence > 0.50) / length(prevalence)
  )

  # Recommendations for PA data
  recommendations <- list(
    analysis_type = "Presence/Absence (PA)",
    primary_function = "fit_longitudinal_PA()",
    rationale = c()
  )

  if (results$diagnostics$prevalence$very_rare > 0.3) {
    recommendations$rationale <- c(recommendations$rationale,
      paste0(round(results$diagnostics$prevalence$very_rare * 100, 1),
             "% of features are very rare (<5% prevalence). PA analysis is appropriate.")
    )
  }

  if (results$diagnostics$prevalence$mean < 0.3) {
    recommendations$rationale <- c(recommendations$rationale,
      paste0("Mean prevalence is ", round(results$diagnostics$prevalence$mean * 100, 1),
             "%. Data is quite sparse - PA analysis recommended.")
    )
  }

  recommendations$alternative <- "If you have the original count data, consider also running abundance analysis with fit_longitudinal_abundance() for complementary insights."

  results$recommendations <- recommendations

  # Make plots
  if (make_plots) {
    # Prevalence distribution
    prevalence_df <- data.frame(prevalence = prevalence)

    p1 <- ggplot2::ggplot(prevalence_df, ggplot2::aes(x = prevalence)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Distribution of Feature Prevalence",
        x = "Prevalence (proportion of samples where feature is present)",
        y = "Number of features"
      ) +
      ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0.10, linetype = "dashed", color = "orange", alpha = 0.5)

    results$plots$prevalence_distribution <- p1
  }

  return(results)
}

#' Diagnose count data
#' @keywords internal
#' @noRd
diagnose_count_data <- function(data_matrix, results, make_plots, n_features_plot) {

  # Calculate key metrics per feature
  feature_means <- rowMeans(data_matrix, na.rm = TRUE)
  feature_vars <- apply(data_matrix, 1, var, na.rm = TRUE)
  zero_props <- rowMeans(data_matrix == 0, na.rm = TRUE)

  # Overall statistics
  overall_mean <- mean(feature_means)
  overall_var <- mean(feature_vars)
  overall_zero_prop <- mean(zero_props)

  # Overdispersion test (variance > mean)
  overdispersion_ratio <- feature_vars / (feature_means + 1e-10)  # Avoid division by zero
  mean_overdispersion <- mean(overdispersion_ratio[is.finite(overdispersion_ratio)])

  results$diagnostics$basic_stats <- list(
    mean_abundance = overall_mean,
    mean_variance = overall_var,
    mean_zero_proportion = overall_zero_prop,
    features_with_high_zeros = sum(zero_props > 0.5) / length(zero_props)
  )

  results$diagnostics$overdispersion <- list(
    mean_ratio = mean_overdispersion,
    prop_overdispersed = sum(overdispersion_ratio > 1, na.rm = TRUE) / sum(is.finite(overdispersion_ratio)),
    interpretation = ifelse(mean_overdispersion > 2, "Strong overdispersion",
                           ifelse(mean_overdispersion > 1.2, "Moderate overdispersion", "Minimal overdispersion"))
  )

  # Zero-inflation test
  # For Poisson, expected zeros = exp(-lambda)
  expected_zeros <- exp(-feature_means)
  excess_zeros <- zero_props - expected_zeros

  results$diagnostics$zero_inflation <- list(
    mean_observed_zeros = overall_zero_prop,
    mean_expected_zeros = mean(expected_zeros, na.rm = TRUE),
    mean_excess_zeros = mean(excess_zeros, na.rm = TRUE),
    prop_excess_zeros = sum(excess_zeros > 0.1, na.rm = TRUE) / length(excess_zeros),
    interpretation = ifelse(mean(excess_zeros, na.rm = TRUE) > 0.15, "Strong zero-inflation",
                           ifelse(mean(excess_zeros, na.rm = TRUE) > 0.05, "Moderate zero-inflation", "Minimal zero-inflation"))
  )

  # Generate recommendations
  recommendations <- list(
    analysis_type = "Count/Abundance",
    primary_function = "fit_longitudinal_abundance()",
    recommended_families = c(),
    rationale = c()
  )

  # Decision logic
  if (overall_zero_prop > 0.7) {
    recommendations$recommended_families <- c("PA analysis (fit_longitudinal_PA)")
    recommendations$rationale <- c(
      paste0("Data is very sparse (", round(overall_zero_prop * 100, 1), "% zeros)."),
      "Consider converting to presence/absence with create_PA_matrix() or fit_longitudinal_PA()."
    )
  }

  if (mean_overdispersion > 1.5) {
    if (results$diagnostics$zero_inflation$interpretation %in% c("Strong zero-inflation", "Moderate zero-inflation")) {
      recommendations$recommended_families <- c(recommendations$recommended_families, "zinbinom2", "zinbinom1")
      recommendations$rationale <- c(recommendations$rationale,
        paste0("Strong overdispersion (ratio: ", round(mean_overdispersion, 2), ") AND zero-inflation detected."),
        "Recommend: family = 'zinbinom2' (zero-inflated negative binomial)"
      )
    } else {
      recommendations$recommended_families <- c(recommendations$recommended_families, "nbinom2", "nbinom1")
      recommendations$rationale <- c(recommendations$rationale,
        paste0("Strong overdispersion detected (ratio: ", round(mean_overdispersion, 2), ")."),
        "Recommend: family = 'nbinom2' (negative binomial, default)"
      )
    }
  } else {
    if (results$diagnostics$zero_inflation$interpretation %in% c("Strong zero-inflation", "Moderate zero-inflation")) {
      recommendations$recommended_families <- c(recommendations$recommended_families, "zipoisson")
      recommendations$rationale <- c(recommendations$rationale,
        "Zero-inflation detected without strong overdispersion.",
        "Recommend: family = 'zipoisson' (zero-inflated Poisson)"
      )
    } else {
      recommendations$recommended_families <- c(recommendations$recommended_families, "poisson", "nbinom2")
      recommendations$rationale <- c(recommendations$rationale,
        "Minimal overdispersion and zero-inflation.",
        "Recommend: family = 'poisson' or 'nbinom2' (try both)"
      )
    }
  }

  recommendations$next_steps <- c(
    "1. Start with the recommended family",
    "2. Check convergence rate (aim for >80%)",
    "3. Compare with alternative families if unsure",
    "4. Consider PA analysis as complementary approach"
  )

  results$recommendations <- recommendations

  # Make plots
  if (make_plots) {
    # 1. Zero proportion distribution
    zero_df <- data.frame(zero_proportion = zero_props)

    p1 <- ggplot2::ggplot(zero_df, ggplot2::aes(x = zero_proportion)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Distribution of Zero Proportions Across Features",
        x = "Proportion of zeros",
        y = "Number of features"
      ) +
      ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5)

    results$plots$zero_distribution <- p1

    # 2. Mean-variance relationship
    mv_df <- data.frame(
      mean = feature_means,
      variance = feature_vars
    ) %>%
      dplyr::filter(mean > 0, variance > 0)  # Log scale needs positive values

    p2 <- ggplot2::ggplot(mv_df, ggplot2::aes(x = mean, y = variance)) +
      ggplot2::geom_point(alpha = 0.3, color = "steelblue") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Mean-Variance Relationship",
        subtitle = "Red line = Poisson (variance = mean). Points above = overdispersion",
        x = "Mean (log scale)",
        y = "Variance (log scale)"
      )

    results$plots$mean_variance <- p2

    # 3. Sample distributions for random features
    set.seed(123)
    n_plot <- min(n_features_plot, nrow(data_matrix))
    sample_features <- sample(nrow(data_matrix), n_plot)

    plot_data <- data_matrix[sample_features, , drop = FALSE] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("feature_id") %>%
      tidyr::pivot_longer(-feature_id, names_to = "sample", values_to = "count") %>%
      dplyr::mutate(feature_id = factor(feature_id, levels = unique(feature_id)))

    p3 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = count)) +
      ggplot2::geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7) +
      ggplot2::facet_wrap(~ feature_id, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Count Distributions for", n_plot, "Random Features"),
        x = "Count",
        y = "Frequency"
      )

    results$plots$example_distributions <- p3
  }

  return(results)
}

#' Print method for longitudinal_diagnosis
#' @export
#' @keywords internal
print.longitudinal_diagnosis <- function(x, ...) {
  cat("==============================================\n")
  cat("LongitudinalSmarter Data Diagnosis\n")
  cat("==============================================\n\n")

  cat("Data Type:", x$data_type, "\n")
  cat("Features:", x$summary_stats$n_features, "\n")
  cat("Samples:", x$summary_stats$n_samples, "\n\n")

  cat("----------------------------------------------\n")
  cat("RECOMMENDATIONS\n")
  cat("----------------------------------------------\n")
  cat("Analysis Type:", x$recommendations$analysis_type, "\n")
  cat("Primary Function:", x$recommendations$primary_function, "\n\n")

  if (length(x$recommendations$recommended_families) > 0) {
    cat("Recommended Model Families:\n")
    for (fam in x$recommendations$recommended_families) {
      cat("  -", fam, "\n")
    }
    cat("\n")
  }

  cat("Rationale:\n")
  for (r in x$recommendations$rationale) {
    cat("  •", r, "\n")
  }

  if (!is.null(x$recommendations$next_steps)) {
    cat("\nNext Steps:\n")
    for (step in x$recommendations$next_steps) {
      cat(" ", step, "\n")
    }
  }

  if (!is.null(x$recommendations$alternative)) {
    cat("\nAlternative:\n")
    cat(" ", x$recommendations$alternative, "\n")
  }

  cat("\n----------------------------------------------\n")
  cat("DIAGNOSTICS SUMMARY\n")
  cat("----------------------------------------------\n")

  if (x$data_type == "PA") {
    cat("Mean Prevalence:", round(x$diagnostics$prevalence$mean * 100, 1), "%\n")
    cat("Very Rare Features (<5%):", round(x$diagnostics$prevalence$very_rare * 100, 1), "%\n")
    cat("Common Features (>50%):", round(x$diagnostics$prevalence$common * 100, 1), "%\n")
  } else {
    cat("Mean Zero Proportion:", round(x$diagnostics$basic_stats$mean_zero_proportion * 100, 1), "%\n")
    cat("Overdispersion:", x$diagnostics$overdispersion$interpretation,
        "(ratio:", round(x$diagnostics$overdispersion$mean_ratio, 2), ")\n")
    cat("Zero-Inflation:", x$diagnostics$zero_inflation$interpretation, "\n")
  }

  cat("\n==============================================\n")
  cat("Use $plots to view diagnostic visualizations\n")
  cat("Use $diagnostics for detailed metrics\n")
  cat("==============================================\n")

  invisible(x)
}
