#' Fit longitudinal presence/absence models with comprehensive contrasts
#'
#' Fits mixed-effects logistic regression models for each feature and extracts:
#' 1. Main model terms (reference group trajectory + interactions)
#' 2. Non-reference group trajectory (vs baseline)
#' 3. Between-group contrasts at each timepoint
#'
#' @param PA_matrix Matrix or data.frame of presence/absence data (features x samples).
#'   Rownames should be feature IDs, colnames should be sample IDs.
#' @param metadata Data frame with sample metadata. Rownames should match colnames of PA_matrix.
#' @param dx_var Name of diagnosis/group variable in metadata (e.g., "Dx.Status")
#' @param timepoint_var Name of timepoint variable in metadata (e.g., "onset_timeline_combined")
#' @param covariates Character vector of covariate names to include in model
#' @param random_var Name of random effect grouping variable (e.g., "patientID")
#' @param baseline_level Baseline timepoint level (default: "t0")
#' @param ref_level Reference level for diagnosis variable (default: "CONTROL")
#' @param non_ref_level Non-reference level for diagnosis variable (default: "CELIAC")
#' @param adjust_method Method for p-value adjustment (default: "BH" for Benjamini-Hochberg)
#' @param return_models Logical; if TRUE, return fitted model objects (warning: memory intensive)
#' @param verbose Logical; if TRUE, print progress messages
#'
#' @return A list with components:
#' \describe{
#'   \item{main_model}{Tibble with main model results (reference group + interactions)}
#'   \item{non_reference}{Tibble with non-reference group trajectory vs baseline}
#'   \item{between_group}{Tibble with between-group contrasts at each timepoint}
#'   \item{combined}{Tibble combining all results with harmonized columns}
#'   \item{n_features}{Number of features analyzed}
#'   \item{n_converged}{Number of models that converged successfully}
#'   \item{failed_features}{Feature IDs where models failed to converge}
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% group_by mutate filter select bind_rows arrange left_join ungroup any_of
#' @importFrom tidyr pivot_longer nest unnest
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom purrr map map_lgl
#' @importFrom glmmTMB glmmTMB
#' @importFrom broom.mixed tidy
#'
#' @examples
#' \dontrun{
#' results <- fit_longitudinal_PA(
#'   PA_matrix = my_PA_data,
#'   metadata = my_metadata,
#'   dx_var = "Dx.Status",
#'   timepoint_var = "onset_timeline_combined",
#'   covariates = c("Country", "Sex", "Age"),
#'   random_var = "patientID"
#' )
#'
#' # View significant results
#' results$combined %>% filter(p.adj < 0.05)
#' }
fit_longitudinal_PA <- function(PA_matrix,
                                metadata,
                                dx_var,
                                timepoint_var,
                                covariates = NULL,
                                random_var = NULL,
                                baseline_level = "t0",
                                ref_level = "CONTROL",
                                non_ref_level = "CELIAC",
                                adjust_method = "BH",
                                return_models = FALSE,
                                verbose = TRUE) {

  # Validate inputs
  validate_inputs(PA_matrix, metadata, dx_var, timepoint_var, covariates, random_var)

  # Remove all-zero features
  feature_sums <- rowSums(PA_matrix)
  PA_matrix <- PA_matrix[feature_sums > 0, , drop = FALSE]

  if (verbose) {
    message("Analyzing ", nrow(PA_matrix), " features across ", ncol(PA_matrix), " samples")
  }

  # Align metadata to PA_matrix columns
  metadata <- metadata[colnames(PA_matrix), , drop = FALSE]

  # Ensure factors have correct levels
  if (!is.factor(metadata[[dx_var]])) {
    metadata[[dx_var]] <- factor(metadata[[dx_var]])
  }
  if (!is.factor(metadata[[timepoint_var]])) {
    metadata[[timepoint_var]] <- factor(metadata[[timepoint_var]])
  }

  # Get timepoint levels
  tp_levels <- levels(metadata[[timepoint_var]])
  if (!(baseline_level %in% tp_levels)) {
    stop("baseline_level '", baseline_level, "' not found in timepoint_var levels")
  }

  # Create long-format data
  if (verbose) message("Converting to long format...")
  pa_long <- PA_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id") %>%
    tidyr::pivot_longer(cols = -feature_id, names_to = "sample_id", values_to = "PA") %>%
    dplyr::left_join(
      metadata %>% tibble::rownames_to_column("sample_id"),
      by = "sample_id"
    )

  # Build formula
  if (is.null(random_var)) {
    formula_str <- paste0("PA ~ ", dx_var, " * ", timepoint_var)
  } else {
    formula_str <- paste0("PA ~ ", dx_var, " * ", timepoint_var, " + (1 | ", random_var, ")")
  }

  if (!is.null(covariates) && length(covariates) > 0) {
    cov_str <- paste(covariates, collapse = " + ")
    formula_str <- gsub("\\+", paste("+", cov_str, "+"), formula_str)
  }

  model_formula <- stats::as.formula(formula_str)

  if (verbose) message("Fitting models with formula: ", formula_str)

  # Fit models for each feature
  if (verbose) message("Fitting models for each feature...")
  fitted_models <- pa_long %>%
    dplyr::group_by(feature_id) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(data, ~ {
        tryCatch({
          glmmTMB::glmmTMB(
            model_formula,
            family = stats::binomial(link = "logit"),
            data = .x
          )
        }, error = function(e) NULL)
      })
    )

  # Count convergence
  n_converged <- sum(!purrr::map_lgl(fitted_models$model, is.null))
  failed_features <- fitted_models$feature_id[purrr::map_lgl(fitted_models$model, is.null)]

  if (verbose) {
    message("Successfully fit ", n_converged, " out of ", nrow(fitted_models), " models")
    if (length(failed_features) > 0) {
      message("Failed to converge: ", length(failed_features), " features")
    }
  }

  # Extract main model results
  if (verbose) message("Extracting main model results...")
  main_results <- fitted_models %>%
    dplyr::mutate(
      tidy_results = purrr::map(model, ~ {
        if (!is.null(.x)) {
          tryCatch(broom.mixed::tidy(.x, conf.int = TRUE), error = function(e) NULL)
        } else NULL
      })
    ) %>%
    dplyr::select(feature_id, tidy_results) %>%
    tidyr::unnest(tidy_results)

  # Apply FDR correction per term
  main_adjusted <- main_results %>%
    dplyr::filter(!is.na(p.value)) %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(p.adj = stats::p.adjust(p.value, method = adjust_method)) %>%
    dplyr::ungroup()

  # Extract non-reference contrasts
  if (verbose) message("Extracting ", non_ref_level, " trajectory contrasts...")
  non_ref_results <- fitted_models %>%
    dplyr::filter(!purrr::map_lgl(model, is.null)) %>%
    dplyr::mutate(
      contrasts = purrr::map(model, ~ extract_non_reference_contrasts(
        .x, tp_levels, dx_var, timepoint_var, non_ref_level, baseline_level
      ))
    ) %>%
    dplyr::select(feature_id, contrasts) %>%
    tidyr::unnest(contrasts)

  # Apply FDR correction per term
  if (nrow(non_ref_results) > 0) {
    non_ref_adjusted <- non_ref_results %>%
      dplyr::filter(!is.na(p.value)) %>%
      dplyr::group_by(term) %>%
      dplyr::mutate(p.adj = stats::p.adjust(p.value, method = adjust_method)) %>%
      dplyr::ungroup()
  } else {
    non_ref_adjusted <- non_ref_results
  }

  # Extract between-group contrasts
  if (verbose) message("Extracting between-group contrasts...")
  between_results <- fitted_models %>%
    dplyr::filter(!purrr::map_lgl(model, is.null)) %>%
    dplyr::mutate(
      contrasts = purrr::map(model, ~ extract_between_group_contrasts(
        .x, dx_var, timepoint_var, ref_level, non_ref_level
      ))
    ) %>%
    dplyr::select(feature_id, contrasts) %>%
    tidyr::unnest(contrasts)

  # Apply FDR correction per timepoint
  if (nrow(between_results) > 0) {
    between_adjusted <- between_results %>%
      dplyr::filter(!is.na(p.value)) %>%
      dplyr::group_by(timepoint) %>%
      dplyr::mutate(p.adj = stats::p.adjust(p.value, method = adjust_method)) %>%
      dplyr::ungroup()
  } else {
    between_adjusted <- between_results
  }

  # Harmonize columns for combined table
  common_spec <- list(
    feature_id = NA_character_,
    term = NA_character_,
    estimate = NA_real_,
    std.error = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_,
    p.adj = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    OR = NA_real_,
    source = NA_character_
  )

  main_std <- main_adjusted %>%
    add_missing_cols_typed(common_spec) %>%
    dplyr::mutate(source = "main_model") %>%
    dplyr::select(names(common_spec))

  non_ref_std <- non_ref_adjusted %>%
    add_missing_cols_typed(common_spec) %>%
    dplyr::mutate(source = "non_reference") %>%
    dplyr::select(names(common_spec))

  between_std <- between_adjusted %>%
    dplyr::select(-dplyr::any_of("timepoint")) %>%
    add_missing_cols_typed(common_spec) %>%
    dplyr::mutate(source = "between_group") %>%
    dplyr::select(names(common_spec))

  # Combine all results
  if (verbose) message("Combining results...")
  combined <- dplyr::bind_rows(main_std, non_ref_std, between_std) %>%
    dplyr::arrange(feature_id, term)

  # Prepare return object
  result <- list(
    main_model = main_adjusted,
    non_reference = non_ref_adjusted,
    between_group = between_adjusted,
    combined = combined,
    n_features = nrow(PA_matrix),
    n_converged = n_converged,
    failed_features = failed_features
  )

  if (return_models) {
    result$models <- fitted_models %>%
      dplyr::filter(!purrr::map_lgl(model, is.null))
  }

  if (verbose) message("Done! ", sum(combined$p.adj < 0.05, na.rm = TRUE), " significant results (FDR < 0.05)")

  return(result)
}
