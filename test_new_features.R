#!/usr/bin/env Rscript

# Test new features in v0.2.0
cat("========================================\n")
cat("Testing LongitudinalSmarter v0.2.0\n")
cat("========================================\n\n")

# Load package
devtools::load_all()

# Create test data
set.seed(123)
cat("Creating test count data...\n")

# 10 features, 20 samples with count data
count_test <- matrix(
  rpois(200, lambda = 5),  # Poisson counts
  nrow = 10,
  ncol = 20,
  dimnames = list(
    paste0("feature_", 1:10),
    paste0("sample_", 1:20)
  )
)

# Metadata
metadata_test <- data.frame(
  Dx.Status = factor(rep(c("CONTROL", "CELIAC"), each = 10)),
  timepoint = factor(rep(c("t0", "t0-6", "t0-12", "t0-18"), 5)),
  patientID = paste0("P", rep(1:10, each = 2)),
  Sex = factor(sample(c("M", "F"), 20, replace = TRUE)),
  row.names = paste0("sample_", 1:20)
)

cat("✓ Test data created\n")
cat("  Features:", nrow(count_test), "\n")
cat("  Samples:", ncol(count_test), "\n")
cat("  Total counts:", sum(count_test), "\n\n")

# Test 1: create_PA_matrix
cat("========================================\n")
cat("Test 1: create_PA_matrix()\n")
cat("========================================\n")
tryCatch({
  PA_from_counts <- create_PA_matrix(
    count_matrix = count_test,
    target_cpm = 0.5,
    min_reads_floor = 3
  )
  cat("✓ create_PA_matrix() works!\n")
  cat("  PA matrix dimensions:", dim(PA_from_counts), "\n")
  cat("  Presence rate:", round(100 * mean(PA_from_counts), 2), "%\n\n")
}, error = function(e) {
  cat("✗ Error in create_PA_matrix():\n")
  cat(conditionMessage(e), "\n\n")
})

# Test 2: fit_longitudinal_abundance with nbinom2
cat("========================================\n")
cat("Test 2: fit_longitudinal_abundance() - nbinom2\n")
cat("========================================\n")
tryCatch({
  results_nb <- fit_longitudinal_abundance(
    count_matrix = count_test,
    metadata = metadata_test,
    dx_var = "Dx.Status",
    timepoint_var = "timepoint",
    covariates = "Sex",
    random_var = "patientID",
    family = "nbinom2",
    verbose = TRUE
  )

  cat("\n✓ Negative binomial model works!\n")
  cat("  Converged models:", results_nb$n_converged, "\n")
  cat("  Results rows:", nrow(results_nb$combined), "\n")
  cat("  Significant (FDR < 0.05):", sum(results_nb$combined$p.adj < 0.05, na.rm = TRUE), "\n\n")
}, error = function(e) {
  cat("\n✗ Error in fit_longitudinal_abundance(nbinom2):\n")
  cat(conditionMessage(e), "\n\n")
})

# Test 3: fit_longitudinal_abundance with Poisson
cat("========================================\n")
cat("Test 3: fit_longitudinal_abundance() - poisson\n")
cat("========================================\n")
tryCatch({
  results_pois <- fit_longitudinal_abundance(
    count_matrix = count_test,
    metadata = metadata_test,
    dx_var = "Dx.Status",
    timepoint_var = "timepoint",
    random_var = "patientID",
    family = "poisson",
    verbose = FALSE
  )

  cat("✓ Poisson model works!\n")
  cat("  Converged models:", results_pois$n_converged, "\n\n")
}, error = function(e) {
  cat("✗ Error in fit_longitudinal_abundance(poisson):\n")
  cat(conditionMessage(e), "\n\n")
})

# Test 4: Zero-inflated Poisson
cat("========================================\n")
cat("Test 4: fit_longitudinal_abundance() - zipoisson\n")
cat("========================================\n")

# Add some zeros to make it more zero-inflated
count_zi <- count_test
count_zi[sample(length(count_zi), 50)] <- 0

tryCatch({
  results_zi <- fit_longitudinal_abundance(
    count_matrix = count_zi,
    metadata = metadata_test,
    dx_var = "Dx.Status",
    timepoint_var = "timepoint",
    random_var = "patientID",
    family = "zipoisson",
    verbose = FALSE
  )

  cat("✓ Zero-inflated Poisson model works!\n")
  cat("  Converged models:", results_zi$n_converged, "\n\n")
}, error = function(e) {
  cat("✗ Error in fit_longitudinal_abundance(zipoisson):\n")
  cat(conditionMessage(e), "\n\n")
})

# Test 5: Check that IRR is calculated
cat("========================================\n")
cat("Test 5: Checking IRR calculation\n")
cat("========================================\n")
if (exists("results_nb")) {
  if ("IRR" %in% names(results_nb$combined)) {
    n_irr <- sum(!is.na(results_nb$combined$IRR))
    cat("✓ IRR column present!\n")
    cat("  Non-NA IRR values:", n_irr, "\n")
    if (n_irr > 0) {
      cat("  Example IRR values:", head(results_nb$combined$IRR[!is.na(results_nb$combined$IRR)], 3), "\n\n")
    }
  } else {
    cat("✗ IRR column missing!\n\n")
  }
}

cat("========================================\n")
cat("All tests complete!\n")
cat("========================================\n")
