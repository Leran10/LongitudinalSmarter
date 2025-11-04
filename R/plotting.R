#' Create heatmaps for longitudinal results by term family
#'
#' Auto-detects different types of contrasts in the results and creates
#' heatmaps showing significant results. Only cells with FDR < threshold are colored.
#'
#' @param df Data frame of results (typically the `combined` output from fit_longitudinal_PA)
#' @param tp_levels Character vector of timepoint levels (excluding baseline if desired)
#' @param fdr FDR threshold for significance (default: 0.05)
#' @param value Type of value to plot: "log_odds" or "OR" (odds ratio)
#' @param min_features Minimum number of significant timepoints required for a feature to be shown
#' @param feature_id_col Name of feature ID column (default: "feature_id")
#' @param ref_level Reference diagnosis level (e.g., "CONTROL")
#' @param non_ref_level Non-reference diagnosis level (e.g., "CELIAC")
#' @param timepoint_var Name used in term labels for timepoints (default: "onset_timeline_combined")
#'
#' @return Named list of ggplot objects, one per detected term family
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 scale_x_discrete scale_y_discrete labs theme_minimal theme element_blank element_text
#' @importFrom dplyr filter mutate group_by summarise arrange pull count %>%
#' @importFrom stringr str_detect str_match
#' @importFrom forcats fct_drop
#'
#' @examples
#' \dontrun{
#' # After fitting models
#' results <- fit_longitudinal_PA(...)
#'
#' # Create heatmaps
#' plots <- plot_longitudinal_heatmaps(
#'   df = results$combined,
#'   tp_levels = c("t0-6", "t0-12", "t0-18", "t0-24"),
#'   fdr = 0.05,
#'   value = "log_odds"
#' )
#'
#' # View individual plots
#' plots$`TP (CONTROL)`
#' plots$`TP (CELIAC)`
#' }
plot_longitudinal_heatmaps <- function(df,
                                       tp_levels,
                                       fdr = 0.05,
                                       value = c("log_odds", "OR"),
                                       min_features = 1,
                                       feature_id_col = "feature_id",
                                       ref_level = "CONTROL",
                                       non_ref_level = "CELIAC",
                                       timepoint_var = "onset_timeline_combined") {

  value <- match.arg(value)

  # Build family catalog based on parameters
  family_catalog <- build_family_catalog(ref_level, non_ref_level, timepoint_var)

  # Detect which families exist
  families <- detect_term_families(df, family_catalog)

  if (length(families) == 0) {
    message("No timepoint-based term families detected.")
    return(list())
  }

  # Create heatmap for each family
  plots <- setNames(vector("list", length(families)), families)
  for (fam in families) {
    plots[[fam]] <- make_term_family_heatmap(
      df = df,
      family = fam,
      family_catalog = family_catalog,
      tp_levels = tp_levels,
      fdr = fdr,
      value = value,
      min_features = min_features,
      feature_id_col = feature_id_col
    )
  }

  return(plots)
}


#' Build catalog of term families
#'
#' @param ref_level Reference diagnosis level
#' @param non_ref_level Non-reference diagnosis level
#' @param timepoint_var Timepoint variable name
#'
#' @return Tibble with family patterns
#' @keywords internal
build_family_catalog <- function(ref_level, non_ref_level, timepoint_var) {
  tibble::tribble(
    ~family, ~match_regex, ~tp_regex, ~note,

    paste0("TP (", ref_level, ")"),
    paste0("^", timepoint_var),
    paste0(timepoint_var, "(t0(?:-(?:\\\\d+|over\\\\d+))?)$"),
    paste0("main TP effects in ", ref_level),

    paste0("TP (", non_ref_level, ")"),
    paste0("^TP \\\\(", non_ref_level, "\\\\):"),
    paste0("TP \\\\(", non_ref_level, "\\\\):\\\\s*([^ ]+)\\\\s+vs\\\\s+t0"),
    paste0(non_ref_level, " vs t0"),

    "Between groups",
    paste0("^", ref_level, " vs ", non_ref_level, " @\\\\s*"),
    paste0("^", ref_level, " vs ", non_ref_level, " @\\\\s*(t0(?:-(?:\\\\d+|over\\\\d+))?)$"),
    "between-group at each timepoint",

    "Interactions",
    paste0(":", timepoint_var),
    paste0(timepoint_var, "(t0(?:-(?:\\\\d+|over\\\\d+))?)$"),
    "interaction terms"
  )
}


#' Detect which term families exist in data
#'
#' @param df Data frame with 'term' column
#' @param family_catalog Catalog of family patterns
#'
#' @return Character vector of detected family names
#' @keywords internal
detect_term_families <- function(df, family_catalog) {
  if (!"term" %in% names(df)) {
    stop("Data frame must have a 'term' column")
  }

  family_catalog %>%
    dplyr::rowwise() %>%
    dplyr::mutate(has = any(grepl(match_regex, df$term))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(has) %>%
    dplyr::pull(family)
}


#' Create heatmap for a single term family
#'
#' @param df Results data frame
#' @param family Family name
#' @param family_catalog Catalog of patterns
#' @param tp_levels Timepoint levels
#' @param fdr FDR threshold
#' @param value "log_odds" or "OR"
#' @param min_features Minimum hits required
#' @param feature_id_col Feature ID column name
#'
#' @return ggplot object
#' @keywords internal
make_term_family_heatmap <- function(df, family, family_catalog, tp_levels,
                                     fdr, value, min_features, feature_id_col) {

  spec <- family_catalog %>% dplyr::filter(family == !!family)
  if (nrow(spec) == 0) {
    stop("Unknown family: ", family)
  }

  # Filter to this family's terms
  sel <- df %>% dplyr::filter(stringr::str_detect(term, spec$match_regex[1]))

  if (nrow(sel) == 0) {
    return(create_empty_plot(paste0("No terms found for family '", family, "'")))
  }

  # Extract timepoint from term
  tp_mat <- stringr::str_match(sel$term, spec$tp_regex[1])
  if (ncol(tp_mat) >= 2) {
    sel$TP <- tp_mat[, 2]
  } else {
    sel$TP <- NA_character_
  }

  # Keep only significant cells
  sel <- sel %>%
    dplyr::mutate(
      TP = factor(TP, levels = tp_levels),
      sig = !is.na(p.adj) & p.adj < fdr,
      value_sig = if (value == "OR") exp(estimate) else estimate
    ) %>%
    dplyr::filter(sig, !is.na(TP))

  if (nrow(sel) == 0) {
    return(create_empty_plot(paste0("No significant results for family '", family, "'")))
  }

  # Filter by minimum hits
  keep_ids <- sel %>%
    dplyr::count(.data[[feature_id_col]], name = "hits") %>%
    dplyr::filter(hits >= min_features) %>%
    dplyr::pull(.data[[feature_id_col]])

  sel <- sel %>% dplyr::filter(.data[[feature_id_col]] %in% keep_ids)

  if (nrow(sel) == 0) {
    return(create_empty_plot(
      paste0("No features pass min_features = ", min_features, " for family '", family, "'")
    ))
  }

  # Order features by max effect size
  ord <- sel %>%
    dplyr::group_by(.data[[feature_id_col]]) %>%
    dplyr::summarise(max_abs = max(abs(value_sig), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(max_abs)) %>%
    dplyr::pull(.data[[feature_id_col]])

  sel[[feature_id_col]] <- factor(sel[[feature_id_col]], levels = ord)
  sel$TP <- forcats::fct_drop(sel$TP)

  # Create plot
  fill_scale <- if (value == "log_odds") {
    ggplot2::scale_fill_gradient2(
      name = "Log-odds",
      low = "#2b8cbe",
      mid = "white",
      high = "#d7301f",
      midpoint = 0
    )
  } else {
    ggplot2::scale_fill_gradient2(
      name = "Odds Ratio",
      low = "#2b8cbe",
      mid = "white",
      high = "#d7301f",
      midpoint = 1,
      trans = "log10"
    )
  }

  ggplot2::ggplot(sel, ggplot2::aes(x = TP, y = .data[[feature_id_col]], fill = value_sig)) +
    ggplot2::geom_tile(width = 0.95, height = 0.95, color = "white", linewidth = 0.2) +
    fill_scale +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::scale_y_discrete(drop = TRUE) +
    ggplot2::labs(
      title = paste0("Heatmap: ", family),
      subtitle = paste0("Significant cells (FDR < ", fdr, "); ", nrow(sel), " cells, ",
                       length(unique(sel[[feature_id_col]])), " features"),
      x = "Timepoint",
      y = "Feature"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      legend.position = "right"
    )
}


#' Create empty plot with message
#'
#' @param message Message to display
#' @return ggplot object
#' @keywords internal
create_empty_plot <- function(message) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = message, size = 5) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1)
}
