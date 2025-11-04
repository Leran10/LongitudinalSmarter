#' Extract non-reference group contrasts from fitted model
#'
#' Uses emmeans to extract trajectory of non-reference group (e.g., CELIAC)
#' comparing each timepoint to baseline within that group.
#'
#' @param mod A fitted glmmTMB model object
#' @param tp_levels Character vector of timepoint levels (must include baseline)
#' @param dx_var Name of diagnosis variable in the model
#' @param timepoint_var Name of timepoint variable in the model
#' @param non_ref_level Level of dx_var to extract (e.g., "CELIAC")
#' @param baseline_level Baseline timepoint (e.g., "t0")
#'
#' @return A tibble with columns: term, estimate, std.error, statistic, p.value, conf.low, conf.high, OR
#' @export
#'
#' @importFrom emmeans emmeans contrast
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter mutate select coalesce
#' @importFrom stringr str_remove
extract_non_reference_contrasts <- function(mod, tp_levels, dx_var = "Dx.Status",
                                            timepoint_var = "onset_timeline_combined",
                                            non_ref_level = "CELIAC",
                                            baseline_level = "t0") {

  if (is.null(mod) || !(baseline_level %in% tp_levels)) {
    return(tibble::tibble())
  }

  # Get emmeans across timepoints within each diagnosis group
  emm <- tryCatch(
    emmeans::emmeans(
      mod,
      stats::as.formula(paste0("~ ", timepoint_var, " | ", dx_var)),
      at = stats::setNames(list(tp_levels), timepoint_var),
      type = "link"  # log-odds scale for logistic models
    ),
    error = function(e) NULL
  )

  if (is.null(emm)) return(tibble::tibble())

  # Contrasts: each TP vs baseline (within each diagnosis group)
  ref_idx <- which(tp_levels == baseline_level)[1]
  con <- tryCatch(
    emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_idx),
    error = function(e) NULL
  )

  if (is.null(con)) return(tibble::tibble())

  # Get summary with inference
  s <- tibble::as_tibble(summary(con, infer = c(TRUE, TRUE)))

  # Unify statistic column (could be z.ratio or t.ratio)
  if (!("statistic" %in% names(s))) {
    zr <- if ("z.ratio" %in% names(s)) s$z.ratio else NULL
    tr <- if ("t.ratio" %in% names(s)) s$t.ratio else NULL
    s$statistic <- dplyr::coalesce(zr, tr, rep(NA_real_, nrow(s)))
  }

  # Ensure CI columns exist
  if (!all(c("lower.CL", "upper.CL") %in% names(s))) {
    ci <- tryCatch(tibble::as_tibble(confint(con)), error = function(e) NULL)
    if (!is.null(ci)) {
      byk <- intersect(names(ci), names(s))
      s <- dplyr::left_join(
        s,
        dplyr::select(ci, dplyr::any_of(c(byk, "lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL"))),
        by = byk
      )
      if (!"lower.CL" %in% names(s) && "asymp.LCL" %in% names(s)) s$lower.CL <- s$asymp.LCL
      if (!"upper.CL" %in% names(s) && "asymp.UCL" %in% names(s)) s$upper.CL <- s$asymp.UCL
    }
  }
  if (!"lower.CL" %in% names(s)) s$lower.CL <- NA_real_
  if (!"upper.CL" %in% names(s)) s$upper.CL <- NA_real_

  # Filter to non-reference group and format
  dx_col <- names(s)[grepl(dx_var, names(s), fixed = TRUE)][1]
  if (is.null(dx_col)) dx_col <- dx_var

  s %>%
    dplyr::filter(.data[[dx_col]] == non_ref_level) %>%
    dplyr::mutate(
      TP = stringr::str_remove(contrast, paste0(" - ", baseline_level, "$")),
      term = paste0("TP (", non_ref_level, "): ", TP, " vs ", baseline_level),
      OR = exp(estimate)
    ) %>%
    dplyr::select(
      term,
      estimate,
      std.error = SE,
      statistic,
      p.value,
      conf.low = lower.CL,
      conf.high = upper.CL,
      OR
    ) %>%
    dplyr::filter(!grepl(paste0("^TP \\(.+\\):\\s*", baseline_level, "\\s+vs\\s+", baseline_level, "$"), term))
}


#' Extract between-group contrasts at each timepoint
#'
#' Uses emmeans to compare diagnosis groups at each timepoint separately.
#'
#' @param mod A fitted glmmTMB model object
#' @param dx_var Name of diagnosis variable in the model
#' @param timepoint_var Name of timepoint variable in the model
#' @param ref_level Reference diagnosis level (e.g., "CONTROL")
#' @param non_ref_level Non-reference diagnosis level (e.g., "CELIAC")
#'
#' @return A tibble with columns: timepoint, term, estimate, std.error, statistic, p.value, conf.low, conf.high, OR
#' @export
#'
#' @importFrom emmeans emmeans contrast
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter mutate select left_join any_of
#' @importFrom stringr str_detect
extract_between_group_contrasts <- function(mod, dx_var = "Dx.Status",
                                            timepoint_var = "onset_timeline_combined",
                                            ref_level = "CONTROL",
                                            non_ref_level = "CELIAC") {

  if (is.null(mod)) return(tibble::tibble())

  # Get emmeans for each group at each timepoint
  emm <- tryCatch(
    emmeans::emmeans(
      mod,
      stats::as.formula(paste0("~ ", dx_var, " | ", timepoint_var))
    ),
    error = function(e) NULL
  )

  if (is.null(emm)) return(tibble::tibble())

  # Contrasts between groups at each timepoint
  contr <- tryCatch(
    emmeans::contrast(emm, method = "revpairwise", by = timepoint_var, adjust = "none"),
    error = function(e) NULL
  )

  if (is.null(contr)) return(tibble::tibble())

  # Get summary and confidence intervals
  s <- tibble::as_tibble(summary(contr))
  ci <- tibble::as_tibble(confint(contr)) %>%
    dplyr::select(dplyr::any_of(c("contrast", timepoint_var, "lower.CL", "upper.CL")))

  keys <- intersect(names(s), names(ci))
  if (length(keys) > 0) {
    s <- dplyr::left_join(s, ci, by = keys)
  }

  # Filter for the contrast we want (non_ref - ref) and flip it to (ref vs non_ref)
  between_df <- s %>%
    dplyr::filter(stringr::str_detect(contrast, paste0("^", non_ref_level, "\\s*-\\s*", ref_level, "$")))

  if (nrow(between_df) == 0) return(tibble::tibble())

  # Check which columns exist
  has_z <- "z.ratio" %in% names(between_df)
  has_t <- "t.ratio" %in% names(between_df)
  has_SE <- "SE" %in% names(between_df)
  has_l <- "lower.CL" %in% names(between_df)
  has_u <- "upper.CL" %in% names(between_df)

  # Get timepoint column name
  tp_col <- names(between_df)[grepl(timepoint_var, names(between_df), fixed = TRUE)][1]
  if (is.null(tp_col)) tp_col <- timepoint_var

  # Flip the contrast: non_ref - ref becomes ref vs non_ref
  between_df %>%
    dplyr::mutate(
      est_flip = -estimate,
      se_use = if (has_SE) SE else NA_real_,
      stat_use = if (has_z) z.ratio else if (has_t) t.ratio else NA_real_,
      low_flip = if (has_l) -upper.CL else NA_real_,
      high_flip = if (has_u) -lower.CL else NA_real_,
      OR = exp(est_flip),
      term = paste0(ref_level, " vs ", non_ref_level, " @ ", .data[[tp_col]]),
      timepoint = .data[[tp_col]]
    ) %>%
    dplyr::select(
      timepoint,
      term,
      estimate = est_flip,
      std.error = se_use,
      statistic = stat_use,
      p.value,
      conf.low = low_flip,
      conf.high = high_flip,
      OR
    )
}
