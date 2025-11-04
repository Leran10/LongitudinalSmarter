# LongitudinalSmarter

> Comprehensive longitudinal mixed-effects models for microbiome data

## Overview

**LongitudinalSmarter** is an R package that streamlines the analysis of longitudinal microbiome data, with a focus on comparing disease vs. control groups over time. The package automatically fits mixed-effects models and extracts multiple types of contrasts that are typically needed but tedious to calculate manually.

### What does it do?

Given presence/absence or abundance data across timepoints and diagnosis groups, the package:

1. **Fits mixed-effects models** for each feature (gene/contig/OTU/ASV)
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

```r
library(LongitudinalSmarter)

# Fit comprehensive longitudinal model
results <- fit_longitudinal_PA(
  PA_matrix = my_PA_data,               # features × samples
  metadata = my_metadata,                # samples × variables
  dx_var = "Dx.Status",                  # "CONTROL" vs "CELIAC"
  timepoint_var = "onset_timeline_combined",  # "t0", "t0-6", etc.
  covariates = c("Country", "Sex", "Age"),
  random_var = "patientID"
)

# View significant results
results$combined %>% filter(p.adj < 0.05)

# Generate heatmaps
plots <- plot_longitudinal_heatmaps(
  df = results$combined,
  tp_levels = c("t0-6", "t0-12", "t0-18", "t0-24"),
  fdr = 0.05
)
```

---

## Features

### Comprehensive Contrast Extraction

The package fits a single mixed-effects model per feature:

```
PA ~ Dx.Status × Timepoint + Covariates + (1|PatientID)
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

The `fit_longitudinal_PA()` function returns a list with:

```r
results$main_model       # Main model terms (reference + interactions)
results$non_reference    # Non-reference group trajectory
results$between_group    # Between-group contrasts
results$combined         # All results combined
results$n_features       # Number of features analyzed
results$n_converged      # Number of successful model fits
results$failed_features  # Features where models failed
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
- `OR`: Odds ratio (for PA models)
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

### Example 4: Custom Heatmaps

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

```
glmmTMB(PA ~ Dx.Status * Timepoint + Covariates + (1|PatientID),
        family = binomial(link = "logit"))
```

- **Family**: Binomial with logit link
- **Fixed effects**: Diagnosis, timepoint, interaction, covariates
- **Random effect**: Patient-specific intercept
- **Estimates**: Log-odds (exponentiate for odds ratios)

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
