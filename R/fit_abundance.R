#' Fit longitudinal abundance models with comprehensive contrasts
#'
#' Fits mixed-effects negative binomial regression models for count data and extracts:
#' 1. Main model terms (reference group trajectory + interactions)
#' 2. Non-reference group trajectory (vs baseline)
#' 3. Between-group contrasts at each timepoint
#'
#' @param count_matrix Matrix or data.frame of count data (features x samples).
#'   Rownames should be feature IDs, colnames should be sample IDs.
#' @param metadata Data frame with sample metadata. Rownames should match colnames of count_matrix.
#' @param dx_var Name of diagnosis/group variable in metadata (e.g., "Dx.Status")
#' @param timepoint_var Name of timepoint variable in metadata (e.g., "onset_timeline_combined")
#' @param covariates Character vector of covariate names to include in model
#' @param random_var Name of random effect grouping variable (e.g., "patientID")
#' @param baseline_level Baseline timepoint level (default: "t0")
#' @param ref_level Reference level for diagnosis variable (default: "CONTROL")
#' @param non_ref_level Non-reference level for diagnosis variable (default: "CELIAC")
#' @param family Model family: "nbinom2" (default), "nbinom1", "poisson",
#'   "truncated_nbinom2", "zipoisson" (zero-inflated Poisson),
#'   "zinbinom2" (zero-inflated negative binomial type 2),
#'   "zinbinom1" (zero-inflated negative binomial type 1)
#' @param adjust_method Method for p-value adjustment (default: "BH")
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
#' @importFrom glmmTMB glmmTMB nbinom2 nbinom1
#' @importFrom broom.mixed tidy
#'
#' @examples
#' \dontrun{
#' results <- fit_longitudinal_abundance(
#'   count_matrix = my_count_data,
#'   metadata = my_metadata,
#'   dx_var = "Dx.Status",
#'   timepoint_var = "onset_timeline_combined",
#'   covariates = c("Country", "Sex", "Age"),
#'   random_var = "patientID",
#'   family = "nbinom2"
#' )
#'
#' # View significant results
#' results$combined %>% filter(p.adj < 0.05)
#' }
fit_longitudinal_abundance <- function(count_matrix,
                                       metadata,
                                       dx_var,
                                       timepoint_var,
                                       covariates = NULL,
                                       random_var = NULL,
                                       baseline_level = "t0",
                                       ref_level = "CONTROL",
                                       non_ref_level = "CELIAC",
                                       family = c("nbinom2", "nbinom1", "poisson", "truncated_nbinom2",
                                                  "zipoisson", "zinbinom2", "zinbinom1"),
                                       adjust_method = "BH",
                                       return_models = FALSE,
                                       verbose = TRUE) {

  family <- match.arg(family)

  # Validate inputs
  validate_inputs(count_matrix, metadata, dx_var, timepoint_var, covariates, random_var)

  # Remove all-zero features
  feature_sums <- rowSums(count_matrix)
  count_matrix <- count_matrix[feature_sums > 0, , drop = FALSE]

  if (verbose) {
    message("Analyzing ", nrow(count_matrix), " features across ", ncol(count_matrix), " samples")
    message("Using ", family, " family")
  }

  # Align metadata to count_matrix columns
  metadata <- metadata[colnames(count_matrix), , drop = FALSE]

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
  count_long <- count_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id") %>%
    tidyr::pivot_longer(cols = -feature_id, names_to = "sample_id", values_to = "count") %>%
    dplyr::left_join(
      metadata %>% tibble::rownames_to_column("sample_id"),
      by = "sample_id"
    )

  # Build formula
  if (is.null(random_var)) {
    formula_str <- paste0("count ~ ", dx_var, " * ", timepoint_var)
  } else {
    formula_str <- paste0("count ~ ", dx_var, " * ", timepoint_var, " + (1 | ", random_var, ")")
  }

  if (!is.null(covariates) && length(covariates) > 0) {
    cov_str <- paste(covariates, collapse = " + ")
    formula_str <- gsub("\\+", paste("+", cov_str, "+"), formula_str)
  }

  model_formula <- stats::as.formula(formula_str)

  # Determine family object
  family_obj <- switch(family,
    "nbinom2" = glmmTMB::nbinom2(),
    "nbinom1" = glmmTMB::nbinom1(),
    "poisson" = stats::poisson(),
    "truncated_nbinom2" = glmmTMB::truncated_nbinom2(),
    "zipoisson" = stats::poisson(),    # Zero-inflation via ziformula
    "zinbinom2" = glmmTMB::nbinom2(),  # Zero-inflation via ziformula
    "zinbinom1" = glmmTMB::nbinom1()   # Zero-inflation via ziformula
  )

  # For zero-inflated models, add zero-inflation formula
  zi_formula <- if (family %in% c("zipoisson", "zinbinom2", "zinbinom1")) {
    stats::as.formula("~1")  # Constant zero-inflation
  } else {
    NULL
  }

  if (verbose) message("Fitting models with formula: ", formula_str)

  # Fit models for each feature
  if (verbose) message("Fitting models for each feature...")
  fitted_models <- count_long %>%
    dplyr::group_by(feature_id) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(data, ~ {
        tryCatch({
          if (!is.null(zi_formula)) {
            glmmTMB::glmmTMB(
              model_formula,
              ziformula = zi_formula,
              family = family_obj,
              data = .x
            )
          } else {
            glmmTMB::glmmTMB(
              model_formula,
              family = family_obj,
              data = .x
            )
          }
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
  # For abundance: use IRR (incidence rate ratio) instead of OR
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
    IRR = NA_real_,
    source = NA_character_
  )

  # Add IRR to main results
  main_adjusted <- main_adjusted %>%
    dplyr::mutate(IRR = exp(estimate))

  # Replace OR with IRR in contrast results
  if (nrow(non_ref_adjusted) > 0 && "OR" %in% names(non_ref_adjusted)) {
    non_ref_adjusted <- non_ref_adjusted %>%
      dplyr::mutate(IRR = OR) %>%
      dplyr::select(-OR)
  }

  if (nrow(between_adjusted) > 0 && "OR" %in% names(between_adjusted)) {
    between_adjusted <- between_adjusted %>%
      dplyr::mutate(IRR = OR) %>%
      dplyr::select(-OR)
  }

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
    n_features = nrow(count_matrix),
    n_converged = n_converged,
    failed_features = failed_features,
    model_family = family
  )

  if (return_models) {
    result$models <- fitted_models %>%
      dplyr::filter(!purrr::map_lgl(model, is.null))
  }

  if (verbose) message("Done! ", sum(combined$p.adj < 0.05, na.rm = TRUE), " significant results (FDR < 0.05)")

  return(result)
}
