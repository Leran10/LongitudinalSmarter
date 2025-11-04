#!/usr/bin/env Rscript

# Test diagnostic dashboard features
cat("========================================\n")
cat("Testing LongitudinalSmarter v0.3.0\n")
cat("Dashboard & Diagnostic Features\n")
cat("========================================\n\n")

# Load package
devtools::load_all()

# Create test data
set.seed(123)
cat("Creating test datasets...\n")

# Test 1: Count data with overdispersion
count_overdispersed <- matrix(
  rnbinom(200, size = 2, mu = 10),  # Negative binomial (overdispersed)
  nrow = 10,
  ncol = 20,
  dimnames = list(
    paste0("feature_", 1:10),
    paste0("sample_", 1:20)
  )
)

# Test 2: Count data with many zeros (zero-inflated)
count_zi <- matrix(
  rpois(200, lambda = 3),
  nrow = 10,
  ncol = 20,
  dimnames = list(
    paste0("feature_", 1:10),
    paste0("sample_", 1:20)
  )
)
# Add excess zeros
zero_mask <- matrix(rbinom(200, 1, 0.4), nrow = 10, ncol = 20)
count_zi[zero_mask == 1] <- 0

# Test 3: Very sparse count data
count_sparse <- matrix(
  rpois(200, lambda = 1),
  nrow = 10,
  ncol = 20,
  dimnames = list(
    paste0("feature_", 1:10),
    paste0("sample_", 1:20)
  )
)
# Make 70% zeros
zero_mask2 <- matrix(rbinom(200, 1, 0.7), nrow = 10, ncol = 20)
count_sparse[zero_mask2 == 1] <- 0

# Test 4: PA data
PA_data <- (count_overdispersed > 5) * 1

cat("✓ Test datasets created\n\n")

# Test 1: Diagnose overdispersed count data
cat("========================================\n")
cat("Test 1: Overdispersed Count Data\n")
cat("========================================\n")
tryCatch({
  diag1 <- diagnose_data(count_overdispersed, data_type = "auto", make_plots = TRUE)
  cat("\n")
  print(diag1)
  cat("\n✓ Overdispersed data diagnosis complete!\n\n")
}, error = function(e) {
  cat("✗ Error:\n")
  cat(conditionMessage(e), "\n\n")
})

# Test 2: Diagnose zero-inflated data
cat("========================================\n")
cat("Test 2: Zero-Inflated Count Data\n")
cat("========================================\n")
tryCatch({
  diag2 <- diagnose_data(count_zi, data_type = "count", make_plots = TRUE)
  cat("\n")
  print(diag2)
  cat("\n✓ Zero-inflated data diagnosis complete!\n\n")
}, error = function(e) {
  cat("✗ Error:\n")
  cat(conditionMessage(e), "\n\n")
})

# Test 3: Diagnose very sparse data
cat("========================================\n")
cat("Test 3: Very Sparse Count Data\n")
cat("========================================\n")
tryCatch({
  diag3 <- diagnose_data(count_sparse, data_type = "auto", make_plots = TRUE)
  cat("\n")
  print(diag3)
  cat("\n✓ Sparse data diagnosis complete!\n\n")
}, error = function(e) {
  cat("✗ Error:\n")
  cat(conditionMessage(e), "\n\n")
})

# Test 4: Diagnose PA data
cat("========================================\n")
cat("Test 4: Presence/Absence Data\n")
cat("========================================\n")
tryCatch({
  diag4 <- diagnose_data(PA_data, data_type = "auto", make_plots = TRUE)
  cat("\n")
  print(diag4)
  cat("\n✓ PA data diagnosis complete!\n\n")
}, error = function(e) {
  cat("✗ Error:\n")
  cat(conditionMessage(e), "\n\n")
})

# Test 5: Check plot generation
cat("========================================\n")
cat("Test 5: Diagnostic Plots\n")
cat("========================================\n")
if (exists("diag1") && !is.null(diag1$plots)) {
  cat("Available plots for overdispersed data:\n")
  cat("  -", paste(names(diag1$plots), collapse = "\n  - "), "\n")
  cat("✓ Plots generated successfully!\n\n")
}

if (exists("diag4") && !is.null(diag4$plots)) {
  cat("Available plots for PA data:\n")
  cat("  -", paste(names(diag4$plots), collapse = "\n  - "), "\n")
  cat("✓ PA plots generated successfully!\n\n")
}

# Test 6: Print diagnostic details
cat("========================================\n")
cat("Test 6: Detailed Diagnostics Access\n")
cat("========================================\n")
if (exists("diag2")) {
  cat("Diagnostic components:\n")
  cat("  - data_type:", diag2$data_type, "\n")
  cat("  - n_features:", diag2$summary_stats$n_features, "\n")
  cat("  - recommendations available:", !is.null(diag2$recommendations), "\n")
  cat("  - diagnostics available:", !is.null(diag2$diagnostics), "\n")
  cat("  - plots available:", !is.null(diag2$plots), "\n")
  cat("✓ All diagnostic components accessible!\n\n")
}

cat("========================================\n")
cat("All dashboard tests complete!\n")
cat("========================================\n\n")

cat("Note: To launch the interactive dashboard, run:\n")
cat("  launch_model_selector()\n\n")

cat("Or with pre-loaded data:\n")
cat("  launch_model_selector(data_matrix = count_overdispersed)\n\n")
