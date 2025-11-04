#' @importFrom dplyr %>%
NULL

#' Validate input data for longitudinal analysis
#'
#' @param PA_matrix Matrix of presence/absence or abundance data (features x samples)
#' @param metadata Data frame of sample metadata
#' @param dx_var Name of diagnosis/group variable
#' @param timepoint_var Name of timepoint variable
#' @param covariates Vector of covariate names
#' @param random_var Name of random effect variable (e.g., "patientID")
#'
#' @return Invisibly returns TRUE if all checks pass, throws error otherwise
#' @keywords internal
validate_inputs <- function(PA_matrix, metadata, dx_var, timepoint_var,
                            covariates = NULL, random_var = NULL) {

  # Check that PA_matrix is a matrix or data.frame
  if (!is.matrix(PA_matrix) && !is.data.frame(PA_matrix)) {
    stop("PA_matrix must be a matrix or data.frame")
  }

  # Check sample alignment
  if (!all(colnames(PA_matrix) %in% rownames(metadata))) {
    stop("Not all sample IDs in PA_matrix columns are found in metadata rownames")
  }

  # Check required variables exist
  required_vars <- c(dx_var, timepoint_var)
  if (!is.null(covariates)) required_vars <- c(required_vars, covariates)
  if (!is.null(random_var)) required_vars <- c(required_vars, random_var)

  missing_vars <- required_vars[!required_vars %in% colnames(metadata)]
  if (length(missing_vars) > 0) {
    stop("Missing required variables in metadata: ", paste(missing_vars, collapse = ", "))
  }

  # Check for factor variables
  if (!is.factor(metadata[[dx_var]])) {
    warning("Converting dx_var '", dx_var, "' to factor")
  }

  if (!is.factor(metadata[[timepoint_var]])) {
    warning("Converting timepoint_var '", timepoint_var, "' to factor")
  }

  # Check for missing values
  if (any(is.na(metadata[, required_vars]))) {
    n_missing <- sum(is.na(metadata[, required_vars]))
    warning("Found ", n_missing, " missing values in metadata variables")
  }

  # Check for all-zero features
  feature_sums <- rowSums(PA_matrix)
  n_zero <- sum(feature_sums == 0)
  if (n_zero > 0) {
    warning("Found ", n_zero, " features with all zeros - these will be removed")
  }

  invisible(TRUE)
}


#' Add missing columns with typed defaults
#'
#' @param df Data frame to modify
#' @param spec Named list of column names and their default values
#'
#' @return Data frame with added columns
#' @keywords internal
add_missing_cols_typed <- function(df, spec) {
  for (nm in names(spec)) {
    if (!nm %in% names(df)) {
      df[[nm]] <- spec[[nm]]
    }
  }
  df
}
