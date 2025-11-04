# LongitudinalSmarter

> Comprehensive longitudinal mixed-effects models for microbiome data

## Overview

**LongitudinalSmarter** is an R package that streamlines the analysis of longitudinal microbiome data, with a focus on comparing disease vs. control groups over time. The package automatically fits mixed-effects models and extracts multiple types of contrasts that are typically needed but tedious to calculate manually.

### What does it do?

Given **presence/absence** or **count/abundance** data across timepoints and diagnosis groups, the package:

1. **Fits mixed-effects models** for each feature (gene/contig/OTU/ASV)
   - Logistic regression for presence/absence data
   - Negative binomial, Poisson, or zero-inflated models for count data
2. **Extracts reference group trajectory** (e.g., CONTROL over time)
3. **Extracts non-reference group trajectory** (e.g., CELIAC over time)
4. **Extracts between-group contrasts** at each timepoint
5. **Applies appropriate FDR correction** for each contrast type
6. **Combines all results** into one comprehensive table
7. **Generates publication-ready heatmaps** showing significant results

### Why use this package?

- **Saves time**: Replaces hundreds of lines of repetitive code
- **Reduces errors**: Consistent workflow across all analyses
- **Complete picture**: Get all relevant contrasts in one analysis
- **Publication-ready**: Auto-generated heatmaps for significant results
- **Flexible**: Works with any longitudinal study design

---

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("Leran10/LongitudinalSmarter")

# Or install locally
devtools::install("path/to/LongitudinalSmarter")
```

---

## Quick Start

### NEW: Interactive Model Selection Dashboard

Not sure which model to use? Launch the interactive dashboard:

```r
library(LongitudinalSmarter)

# Launch dashboard to upload and diagnose your data
launch_model_selector()

# Or diagnose data programmatically
diagnosis <- diagnose_data(my_count_matrix)
diagnosis  # View recommendations
diagnosis$plots$mean_variance  # View diagnostic plots
```

The dashboard will:
- Analyze your data characteristics (sparsity, overdispersion, zero-inflation)
- Recommend appropriate models with clear rationale
- Generate diagnostic plots
- Provide ready-to-use R code

### Presence/Absence Analysis

```r
library(LongitudinalSmarter)

# Fit comprehensive longitudinal PA model
results_PA <- fit_longitudinal_PA(
  PA_matrix = my_PA_data,               # features × samples (0/1)
  metadata = my_metadata,                # samples × variables
  dx_var = "Dx.Status",                  # "CONTROL" vs "CELIAC"
  timepoint_var = "onset_timeline_combined",  # "t0", "t0-6", etc.
  covariates = c("Country", "Sex", "Age"),
  random_var = "patientID"
)

# View significant results
results_PA$combined %>% filter(p.adj < 0.05)

# Generate heatmaps
plots <- plot_longitudinal_heatmaps(
  df = results_PA$combined,
  tp_levels = c("t0-6", "t0-12", "t0-18", "t0-24"),
  fdr = 0.05
)
```

### Abundance/Count Analysis

```r
# Option 1: Create PA matrix from count data first
PA_matrix <- create_PA_matrix(
  count_matrix = my_count_data,
  target_cpm = 0.5,         # CPM threshold
  min_reads_floor = 3       # Minimum read count
)

# Option 2: Analyze abundance directly with negative binomial models
results_abundance <- fit_longitudinal_abundance(
  count_matrix = my_count_data,         # features × samples (counts)
  metadata = my_metadata,
  dx_var = "Dx.Status",
  timepoint_var = "onset_timeline_combined",
  covariates = c("Country", "Sex", "Age"),
  random_var = "patientID",
  family = "nbinom2"          # negative binomial type 2
)

# Results include IRR (Incidence Rate Ratios) instead of OR
results_abundance$combined %>% filter(p.adj < 0.05)
```

---

## Features

### Two Analysis Types

**1. Presence/Absence (PA) Analysis** - `fit_longitudinal_PA()`
- Binary logistic regression
- Returns Odds Ratios (OR)
- Best for: sparse features, detection/occurrence questions

**2. Abundance/Count Analysis** - `fit_longitudinal_abundance()`
- Negative binomial, Poisson, or zero-inflated models
- Returns Incidence Rate Ratios (IRR)
- Best for: count data, quantitative abundance questions
- Supported families:
  - `nbinom2`: Negative binomial type 2 (default, recommended)
  - `nbinom1`: Negative binomial type 1
  - `poisson`: Poisson regression
  - `truncated_nbinom2`: Truncated negative binomial
  - `zipoisson`: Zero-inflated Poisson
  - `zinbinom2`: Zero-inflated negative binomial type 2
  - `zinbinom1`: Zero-inflated negative binomial type 1

### Diagnostic Tools (NEW in v0.3.0)

**1. Interactive Dashboard** - `launch_model_selector()`
- Upload your data through a web interface
- Get instant model recommendations
- View diagnostic plots interactively
- Copy ready-to-use R code

```r
# Launch the dashboard
launch_model_selector()

# Or pre-load your data
launch_model_selector(data_matrix = my_counts, metadata = my_metadata)
```

**2. Programmatic Diagnostics** - `diagnose_data()`
- Analyze data characteristics automatically
- Get model recommendations with rationale
- Generate diagnostic plots
- Access detailed metrics

```r
# Diagnose your data
diagnosis <- diagnose_data(my_count_matrix)

# View recommendations
print(diagnosis)

# Access diagnostic plots
diagnosis$plots$mean_variance      # Mean-variance relationship
diagnosis$plots$zero_distribution  # Zero proportion distribution
diagnosis$plots$example_distributions  # Sample feature distributions

# Access detailed metrics
diagnosis$diagnostics$overdispersion
diagnosis$diagnostics$zero_inflation
```

**What the diagnostics detect:**
- **Data type**: Automatically identifies PA vs count data
- **Sparsity**: Calculates zero proportions across features
- **Overdispersion**: Tests if variance > mean (nbinom2 vs poisson)
- **Zero-inflation**: Detects excess zeros beyond expected (zero-inflated models)
- **Prevalence**: For PA data, identifies rare vs common features

**Example output:**
```
==============================================
LongitudinalSmarter Data Diagnosis
==============================================
Data Type: count
Features: 1000
Samples: 50

RECOMMENDATIONS
Analysis Type: Count/Abundance
Primary Function: fit_longitudinal_abundance()
Recommended Model Families: nbinom2

Rationale:
  • Strong overdispersion detected (ratio: 5.2)
  • Recommend: family = 'nbinom2' (negative binomial)

Next Steps:
  1. Start with the recommended family
  2. Check convergence rate (aim for >80%)
  3. Compare with alternative families if unsure
==============================================
```

---

## Model Selection Guide

### Should I use PA or Abundance analysis?

**Use Presence/Absence (PA) analysis when:**
- Your data is very sparse (many zeros)
- You care about detection/occurrence rather than quantity
- You want to answer: "Is this feature more/less likely to be present?"
- Features have highly variable counts that don't correlate with biological signal
- You're analyzing amplicon data (16S, ITS) where counts may be unreliable

**Use Abundance analysis when:**
- You have count data with meaningful quantitative information
- You care about changes in abundance, not just presence
- You want to answer: "Does this feature increase/decrease in abundance?"
- You're analyzing shotgun metagenomic or RNA-seq data
- Features show consistent abundance patterns across replicates

**Not sure? Do both!** PA and abundance analyses often reveal complementary insights.

### Which model family for abundance analysis?

#### Quick Decision Tree:

```
START: Do you have count data?
  │
  ├─ YES → Are there many zeros (>40% of values)?
  │   │
  │   ├─ YES → Are zeros more than expected by chance?
  │   │   │
  │   │   ├─ YES → Use zero-inflated models (zipoisson or zinbinom2)
  │   │   └─ NO  → Use nbinom2 (default)
  │   │
  │   └─ NO → Use nbinom2 (default)
  │
  └─ NO → Use fit_longitudinal_PA() instead
```

#### Detailed Model Selection:

**1. Start with nbinom2 (default, recommended for most cases)**
- **When**: Most microbiome count data
- **Why**: Handles overdispersion (variance > mean), which is nearly universal in microbiome data
- **Interpretation**: IRR = multiplicative change in expected count

```r
results <- fit_longitudinal_abundance(
  count_matrix = my_counts,
  metadata = my_metadata,
  dx_var = "Dx.Status",
  timepoint_var = "timepoint",
  random_var = "patientID",
  family = "nbinom2"  # Start here!
)
```

**2. Consider poisson only if:**
- Your mean ≈ variance (very rare in microbiome data)
- You have theoretical reasons to expect equidispersion
- **Check**: Compare nbinom2 vs poisson - if results are similar, data may fit Poisson

**3. Use zero-inflated models (zipoisson, zinbinom2) when:**
- You have **excess zeros** beyond what the count distribution predicts
- Biological zeros (true absence) mixed with sampling zeros (not detected)
- Examples: Host-associated microbes (truly absent in some individuals)

**How to detect excess zeros:**
```r
# Calculate proportion of zeros
zero_prop <- rowMeans(my_counts == 0)
hist(zero_prop, main = "Distribution of zero proportions across features")

# If many features have >50% zeros, consider zero-inflated models
```

**4. Use zinbinom2 over zipoisson when:**
- You have excess zeros AND overdispersion
- Most common choice for zero-inflated microbiome data
- More flexible than zipoisson

**5. Use truncated_nbinom2 when:**
- Zero counts are structurally impossible (e.g., analysis restricted to features that must be present)
- Very rare in practice

**6. nbinom1 vs nbinom2:**
- **nbinom2** (default): Variance = μ + μ²/θ (quadratic mean-variance)
- **nbinom1**: Variance = μ × φ (linear mean-variance)
- Use nbinom1 if you have theoretical reasons, but nbinom2 is more common

### How to compare models empirically

```r
# Fit multiple models
families <- c("nbinom2", "poisson", "zinbinom2")

model_comparison <- lapply(families, function(fam) {
  fit_longitudinal_abundance(
    count_matrix = my_counts,
    metadata = my_metadata,
    dx_var = "Dx.Status",
    timepoint_var = "timepoint",
    random_var = "patientID",
    family = fam,
    verbose = FALSE
  )
})
names(model_comparison) <- families

# Compare convergence (higher is better)
sapply(model_comparison, function(x) x$n_converged)

# Compare number of significant results
sapply(model_comparison, function(x) {
  sum(x$combined$p.adj < 0.05, na.rm = TRUE)
})

# Check for features that are significant in one model but not another
sig_nbinom2 <- model_comparison$nbinom2$combined %>%
  filter(p.adj < 0.05) %>%
  pull(feature_id) %>%
  unique()

sig_zinbinom2 <- model_comparison$zinbinom2$combined %>%
  filter(p.adj < 0.05) %>%
  pull(feature_id) %>%
  unique()

# Features only significant in zero-inflated model
setdiff(sig_zinbinom2, sig_nbinom2)
```

### General Recommendations

**For most microbiome studies:**
1. **Default**: Use `family = "nbinom2"` for abundance
2. **High zeros**: If >50% zeros across many features, try `zinbinom2`
3. **Always check**: Compare convergence rates and inspect failed features

**Red flags:**
- Many convergence failures (>20%) → Data may be too sparse, consider PA analysis
- Huge differences between models → Data may not fit assumptions, interpret cautiously
- All p-values identical across models → Problem with data or model specification

### Reporting your choice

When publishing, report:
1. Which model family you used
2. Why you chose it (e.g., "nbinom2 selected due to overdispersion")
3. Convergence rate (e.g., "95% of features converged")
4. Sensitivity analysis if relevant (e.g., "Results robust to choice of poisson vs nbinom2")

---

### Comprehensive Contrast Extraction

The package fits a single mixed-effects model per feature:

```
PA ~ Dx.Status × Timepoint + Covariates + (1|PatientID)
# or
Count ~ Dx.Status × Timepoint + Covariates + (1|PatientID)
```

From this model, it automatically extracts:

| Contrast Type | Example Terms | Scientific Question |
|--------------|---------------|---------------------|
| **Reference trajectory** | `onset_timeline_combined_t0-6` | How does CONTROL group change over time? |
| **Non-reference trajectory** | `TP (CELIAC): t0-6 vs t0` | How does CELIAC group change over time? |
| **Between-group at each TP** | `CONTROL vs CELIAC @ t0-6` | Do groups differ at this specific timepoint? |
| **Interactions** | `Dx.Status:onset_timeline_combined` | Does the trajectory differ between groups? |

### Smart FDR Correction

Different contrasts are corrected using different strategies:

- **Main model terms**: FDR correction grouped by term type
- **Non-reference contrasts**: FDR correction per timepoint comparison
- **Between-group contrasts**: FDR correction per timepoint

### Auto-Generated Heatmaps

The `plot_longitudinal_heatmaps()` function:

- Auto-detects which types of contrasts exist in your results
- Creates separate heatmaps for each contrast family
- Shows only FDR-significant cells
- Orders features by effect size
- Customizable colors, FDR thresholds, and display options

---

## Input Data Format

### PA Matrix (or Abundance Matrix)

- **Rows**: Features (genes, contigs, OTUs, ASVs, etc.)
- **Columns**: Samples
- **Rownames**: Feature IDs
- **Colnames**: Sample IDs matching metadata rownames
- **Values**: 0/1 for presence/absence, or counts for abundance

```r
# Example PA matrix
#                    Sample1 Sample2 Sample3 Sample4
# feature_1               1       0       1       1
# feature_2               0       0       1       0
# feature_3               1       1       1       1
```

### Metadata

- **Rownames**: Sample IDs matching PA matrix colnames
- **Required columns**:
  - Diagnosis/group variable (e.g., "Dx.Status")
  - Timepoint variable (e.g., "onset_timeline_combined")
  - Random effect variable (e.g., "patientID")
  - Any covariates you want to adjust for

```r
# Example metadata
#         Dx.Status onset_timeline_combined patientID Country Sex
# Sample1 CONTROL   t0                      P001      US      F
# Sample2 CELIAC    t0-6                    P002      US      M
# Sample3 CONTROL   t0-12                   P001      US      F
# Sample4 CELIAC    t0-6                    P003      Italy   F
```

---

## Output Structure

Both `fit_longitudinal_PA()` and `fit_longitudinal_abundance()` return a list with:

```r
results$main_model       # Main model terms (reference + interactions)
results$non_reference    # Non-reference group trajectory
results$between_group    # Between-group contrasts
results$combined         # All results combined
results$n_features       # Number of features analyzed
results$n_converged      # Number of successful model fits
results$failed_features  # Features where models failed
results$model_family     # Model family used
```

Each results table contains:

- `feature_id`: Feature identifier
- `term`: Model term or contrast name
- `estimate`: Log-odds (PA) or log-rate (abundance)
- `std.error`: Standard error
- `statistic`: Test statistic (z or t)
- `p.value`: Raw p-value
- `p.adj`: FDR-adjusted p-value
- `conf.low`, `conf.high`: Confidence interval
- `OR`: Odds ratio (PA models only)
- `IRR`: Incidence rate ratio (abundance models only)
- `source`: Which contrast type

---

## Examples

### Example 1: Basic PA Analysis

```r
# Fit model
results <- fit_longitudinal_PA(
  PA_matrix = total.PA,
  metadata = total.metadata,
  dx_var = "Dx.Status",
  timepoint_var = "onset_timeline_combined",
  covariates = c("Sex", "Age"),
  random_var = "patientID"
)

# How many features show significant changes in CELIAC group?
results$non_reference %>%
  filter(p.adj < 0.05) %>%
  pull(feature_id) %>%
  unique() %>%
  length()
```

### Example 2: Timepoint-Specific Differences

```r
# Which features differ between groups at t0-12?
results$between_group %>%
  filter(timepoint == "t0-12", p.adj < 0.05) %>%
  arrange(p.adj)
```

### Example 3: Strongest Effects

```r
# Top 10 features with largest effect sizes
results$combined %>%
  filter(p.adj < 0.05) %>%
  group_by(feature_id) %>%
  summarise(max_effect = max(abs(estimate))) %>%
  arrange(desc(max_effect)) %>%
  head(10)
```

### Example 4: Abundance Analysis with Different Models

```r
# Compare different model families
models_to_test <- c("nbinom2", "poisson", "zipoisson")

abundance_results <- lapply(models_to_test, function(fam) {
  fit_longitudinal_abundance(
    count_matrix = my_counts,
    metadata = my_metadata,
    dx_var = "Dx.Status",
    timepoint_var = "onset_timeline_combined",
    random_var = "patientID",
    family = fam,
    verbose = FALSE
  )
})
names(abundance_results) <- models_to_test

# Compare convergence rates
sapply(abundance_results, function(x) x$n_converged)

# Compare number of significant results
sapply(abundance_results, function(x) {
  sum(x$combined$p.adj < 0.05, na.rm = TRUE)
})
```

### Example 5: Custom Heatmaps

```r
# Create heatmaps with custom settings
plots <- plot_longitudinal_heatmaps(
  df = results$combined,
  tp_levels = c("t0-6", "t0-12", "t0-18", "t0-24"),
  fdr = 0.01,              # Stricter threshold
  value = "OR",            # Show odds ratios instead of log-odds
  min_features = 2         # Require ≥2 significant timepoints per feature
)

# Save plots
ggsave("heatmap_control.pdf", plots$`TP (CONTROL)`, width = 10, height = 8)
ggsave("heatmap_celiac.pdf", plots$`TP (CELIAC)`, width = 10, height = 8)
```

---

## Model Details

### Mixed-Effects Logistic Regression (PA)

For presence/absence data, the package fits:

```r
glmmTMB(PA ~ Dx.Status * Timepoint + Covariates + (1|PatientID),
        family = binomial(link = "logit"))
```

- **Family**: Binomial with logit link
- **Fixed effects**: Diagnosis, timepoint, interaction, covariates
- **Random effect**: Patient-specific intercept
- **Estimates**: Log-odds (exponentiate for odds ratios)

### Mixed-Effects Models for Count Data (Abundance)

For count/abundance data, the package fits:

```r
glmmTMB(Count ~ Dx.Status * Timepoint + Covariates + (1|PatientID),
        family = nbinom2())  # or other families
```

**Recommended families by data type**:
- **nbinom2** (default): Best for most microbiome count data with overdispersion
- **poisson**: Use when variance ≈ mean (rare in microbiome data)
- **zipoisson/zinbinom2**: Use when excess zeros beyond what's expected from the count distribution
- **truncated_nbinom2**: Use when zero counts are structurally impossible

**Estimates**: Log-rate (exponentiate for incidence rate ratios)

### Helper Function: create_PA_matrix()

Converts count data to presence/absence using TMM normalization:

```r
PA_matrix <- create_PA_matrix(
  count_matrix = my_counts,
  target_cpm = 0.5,           # CPM threshold
  min_reads_floor = 3,        # Minimum raw count
  normalization_method = "TMMwsp"  # TMM for sparse data
)
```

A feature is "present" if **both** conditions are met:
- Normalized CPM ≥ target_cpm
- Raw count ≥ min_reads_floor

This prevents calling rare features "present" due to normalization artifacts.

### Why emmeans?

The package uses **emmeans** to extract contrasts not directly available from model coefficients:

- Non-reference group trajectory (e.g., CELIAC changes vs baseline)
- Between-group differences at each timepoint

This avoids manual re-leveling and re-fitting models.

---

## Troubleshooting

### "Model failed to converge"

Some features may have convergence issues due to:

- Complete separation (feature present in ALL cases or ALL controls)
- Sparse data (feature present in <5% of samples)
- Multicollinearity among covariates

**Solutions**:
- Check `results$failed_features` to see which features failed
- Consider filtering very rare features before analysis
- Check for confounded variables in metadata

### "No significant results"

If FDR-adjusted p-values are all >0.05:

- Try less stringent threshold (e.g., `fdr = 0.1`)
- Check if raw p-values show any signal
- Verify sample size is adequate for your effect sizes
- Ensure metadata variables are correctly formatted

### "Memory issues with large datasets"

For >10,000 features:

- Set `return_models = FALSE` (default) to avoid storing model objects
- Process cohorts separately (US, Italy) then combine results
- Consider filtering features by prevalence/abundance first

---

## Comparison with Original Workflow

| Original Script | LongitudinalSmarter Package |
|----------------|----------------------------|
| ~7000 lines | ~50 lines |
| Manual contrast extraction | Automatic |
| Separate chunks for each cohort | Single function call |
| Manual FDR correction | Automatic, appropriate strategy |
| Custom heatmap code | Auto-generated |
| Hard to modify | Easily customizable |

---

## Citation

If you use this package, please cite:

```
[Citation to be added]
```

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

---

## License

MIT License - see LICENSE file

---

## Contact

- **Issues**: [GitHub Issues](https://github.com/Leran10/LongitudinalSmarter/issues)
- **Author**: Leran Wang
- **Email**: leranwang10@gmail.com

---

## See Also

- **glmmTMB**: [Mixed models package](https://github.com/glmmTMB/glmmTMB)
- **emmeans**: [Estimated marginal means](https://github.com/rvlenth/emmeans)
- **edgeR**: [RNA-seq analysis](https://bioconductor.org/packages/edgeR/)
