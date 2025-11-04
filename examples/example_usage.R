# Example usage of LongitudinalSmarter package
# This script demonstrates how to use the package for analyzing
# longitudinal presence/absence data in celiac disease study

library(LongitudinalSmarter)
library(dplyr)

# =============================================================================
# 1. LOAD YOUR DATA
# =============================================================================

# Load PA matrix (features x samples)
# Rownames = feature IDs, colnames = sample IDs
total.PA <- readRDS("path/to/your/PA_matrix.rds")

# Load metadata (samples as rownames)
total.metadata <- readRDS("path/to/your/metadata.rds")

# =============================================================================
# 2. FIT COMPREHENSIVE LONGITUDINAL MODEL
# =============================================================================

results <- fit_longitudinal_PA(
  PA_matrix = total.PA,
  metadata = total.metadata,
  dx_var = "Dx.Status",                    # Diagnosis variable
  timepoint_var = "onset_timeline_combined", # Timepoint variable
  covariates = c(
    "Country",
    "Sex",
    "Age.at.Gluten.Introduction..months.",
    "HLA.Category",
    "feeding_first_year",
    "Delivery.Mode"
  ),
  random_var = "patientID",                 # Random effect
  baseline_level = "t0",                    # Baseline timepoint
  ref_level = "CONTROL",                    # Reference group
  non_ref_level = "CELIAC",                 # Non-reference group
  adjust_method = "BH",                     # FDR correction method
  return_models = FALSE,                    # Set TRUE to keep model objects
  verbose = TRUE                            # Print progress
)

# =============================================================================
# 3. EXAMINE RESULTS
# =============================================================================

# Check summary
cat("Total features analyzed:", results$n_features, "\n")
cat("Successfully converged:", results$n_converged, "\n")
cat("Failed features:", length(results$failed_features), "\n")

# View structure
str(results, max.level = 1)

# Results components:
# - results$main_model: Reference group trajectory + interactions
# - results$non_reference: Non-reference group (CELIAC) trajectory
# - results$between_group: Group differences at each timepoint
# - results$combined: All results combined

# =============================================================================
# 4. FILTER SIGNIFICANT RESULTS
# =============================================================================

# Overall significant results (FDR < 0.05)
sig_all <- results$combined %>%
  filter(p.adj < 0.05)

cat("\nTotal significant results:", nrow(sig_all), "\n")

# Significant by source
sig_all %>%
  count(source) %>%
  print()

# Significant CELIAC trajectory changes
sig_celiac <- results$non_reference %>%
  filter(p.adj < 0.05)

cat("\nSignificant CELIAC trajectory changes:", nrow(sig_celiac), "\n")
print(head(sig_celiac))

# Significant between-group differences
sig_between <- results$between_group %>%
  filter(p.adj < 0.05)

cat("\nSignificant between-group differences:", nrow(sig_between), "\n")
print(head(sig_between))

# =============================================================================
# 5. CREATE HEATMAPS
# =============================================================================

# Define timepoints (excluding baseline for visualization)
tp_levels <- c("t0-6", "t0-12", "t0-18", "t0-24", "t0-30", "t0-36", "t0-over42")

# Generate heatmaps for all detected term families
plots <- plot_longitudinal_heatmaps(
  df = results$combined,
  tp_levels = tp_levels,
  fdr = 0.05,
  value = "log_odds",     # Can also use "OR" for odds ratios
  min_features = 1,
  feature_id_col = "feature_id",
  ref_level = "CONTROL",
  non_ref_level = "CELIAC",
  timepoint_var = "onset_timeline_combined"
)

# View available plots
names(plots)

# Display individual plots
if ("TP (CONTROL)" %in% names(plots)) {
  print(plots$`TP (CONTROL)`)
}

if ("TP (CELIAC)" %in% names(plots)) {
  print(plots$`TP (CELIAC)`)
}

if ("Between groups" %in% names(plots)) {
  print(plots$`Between groups`)
}

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

# Save combined results table
write.csv(
  results$combined,
  "longitudinal_PA_results_combined.csv",
  row.names = FALSE
)

# Save significant results only
write.csv(
  sig_all,
  "longitudinal_PA_results_significant.csv",
  row.names = FALSE
)

# Save plots
if (length(plots) > 0) {
  for (plot_name in names(plots)) {
    filename <- paste0("heatmap_", gsub(" ", "_", plot_name), ".pdf")
    ggsave(filename, plots[[plot_name]], width = 10, height = 8)
  }
}

# =============================================================================
# 7. ADVANCED: CUSTOM FILTERING AND VISUALIZATION
# =============================================================================

# Find features with significant changes in CELIAC at specific timepoints
celiac_t0_6 <- results$non_reference %>%
  filter(grepl("t0-6", term), p.adj < 0.05) %>%
  arrange(p.adj)

cat("\nFeatures with significant CELIAC changes at t0-6:", nrow(celiac_t0_6), "\n")

# Find features that differ between groups at baseline
baseline_diff <- results$between_group %>%
  filter(timepoint == "t0", p.adj < 0.05) %>%
  arrange(desc(abs(estimate)))

cat("\nFeatures different between groups at baseline:", nrow(baseline_diff), "\n")

# Identify features with strongest effects
top_features <- results$combined %>%
  filter(p.adj < 0.05) %>%
  group_by(feature_id) %>%
  summarise(
    max_effect = max(abs(estimate), na.rm = TRUE),
    n_sig_terms = n()
  ) %>%
  arrange(desc(max_effect)) %>%
  head(20)

cat("\nTop 20 features by effect size:\n")
print(top_features)

# =============================================================================
# 8. EXPORT FOR DOWNSTREAM ANALYSIS
# =============================================================================

# Create wide-format summary for each feature
feature_summary <- results$combined %>%
  group_by(feature_id) %>%
  summarise(
    n_sig_results = sum(p.adj < 0.05, na.rm = TRUE),
    min_pval = min(p.value, na.rm = TRUE),
    min_padj = min(p.adj, na.rm = TRUE),
    max_effect = max(abs(estimate), na.rm = TRUE)
  ) %>%
  filter(n_sig_results > 0) %>%
  arrange(desc(n_sig_results), min_padj)

write.csv(feature_summary, "feature_summary.csv", row.names = FALSE)

cat("\n\nAnalysis complete!\n")
