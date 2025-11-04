# LongitudinalSmarter R Package Development

## Project Overview

**Goal**: Create an R package to streamline longitudinal mixed-effects modeling for microbiome data (phage contigs, bacteria, etc.), particularly for comparing disease vs. control groups over time.

**Original Problem**: The existing analysis script (`Celiac_phage_orf_contig_phrog_PA_abundance_analysis.Rmd`) contains ~7000 lines of repetitive code for analyzing celiac disease vs. control phage data across multiple timepoints.

**Solution**: Package the workflow into reusable functions that:
1. Fit mixed-effects models (presence/absence or abundance)
2. Extract reference group trajectory (e.g., CONTROL over time)
3. Extract non-reference group trajectory (e.g., CELIAC over time)
4. Extract between-group contrasts at each timepoint
5. Combine all results into one comprehensive table
6. Generate publication-ready heatmaps

---

## Current Progress

### ✅ Completed Tasks

#### 1. Package Structure Created
```
LongitudinalSmarter/
├── DESCRIPTION          # Package metadata and dependencies
├── NAMESPACE           # Exported functions
├── LICENSE             # MIT license
├── R/                  # Function source files
│   └── utils.R         # Validation and helper functions
├── man/                # Documentation (to be generated)
├── tests/testthat/     # Unit tests (to be added)
└── data-raw/           # Example data scripts
```

#### 2. DESCRIPTION File
- Package name: `LongitudinalSmarter`
- Dependencies specified:
  - Core: dplyr, tidyr, purrr, tibble, stringr
  - Modeling: glmmTMB, emmeans, broom.mixed, edgeR
  - Visualization: ggplot2, forcats
- Ready for CRAN-style package building

#### 3. Utility Functions (`R/utils.R`)
- `validate_inputs()`: Comprehensive input validation
  - Checks sample alignment between data matrix and metadata
  - Verifies required variables exist
  - Warns about factor conversion needs
  - Detects missing values and all-zero features
- `add_missing_cols_typed()`: Safely adds missing columns with correct types

---

## Original Workflow Analysis

### Chunks Analyzed from Original Script

**Chunk 1: Main Model** (lines 656-701)
```r
# Fits glmmTMB logistic regression per contig
PA ~ Dx.Status * onset_timeline_combined + covariates + (1|patientID)
```
- **Output**: Reference group (CONTROL) trajectory + interaction terms
- **Limitation**: Non-reference group (CELIAC) trajectory not directly available

**Chunk 2: Non-Reference Stats** (lines 704-831)
```r
# Uses emmeans to extract CELIAC trajectory
extract_celiac_tp_vs_t0()
```
- **Output**: CELIAC timepoint contrasts vs. baseline (t0)
- **Terms**: "TP (CELIAC): t0-6 vs t0", "TP (CELIAC): t0-12 vs t0", etc.

**Chunk 3: Between-Group Contrasts** (lines 834-1011)
```r
# Uses emmeans for CONTROL vs CELIAC at each timepoint
extract_contrasts() with method = "revpairwise"
```
- **Output**: Between-group differences at each timepoint
- **Terms**: "CONTROL vs CELIAC @ t0", "CONTROL vs CELIAC @ t0-6", etc.
- **Note**: Flips "CELIAC - CONTROL" to "CONTROL vs CELIAC" for consistency

**Chunk 4: Heatmaps** (lines 1014-1033)
```r
# Auto-detects term families and creates heatmaps
auto_heatmaps_detected()
```
- **Output**: ggplot2 heatmaps showing significant results
- **Features**: Only shows FDR-significant cells, orders by effect size

---

## Package Design Decisions

### 1. Model Extensibility
**Current Scope** (v0.1.0):
- Binary diagnosis groups (e.g., CONTROL vs CELIAC)
- Categorical timepoints with baseline (t0)
- glmmTMB logistic (PA) and negative binomial (abundance)

**Future Extensions** (planned):
- Zero-inflated models
- 3+ diagnosis groups
- Continuous time variables
- Different model frameworks (lme4, brms)

### 2. Input Requirements
**Required inputs**:
- `PA_matrix`: Features × samples matrix (rownames = feature IDs, colnames = sample IDs)
- `metadata`: Data frame with samples as rownames
- `dx_var`: Column name for diagnosis/group variable
- `timepoint_var`: Column name for timepoint variable
- `covariates`: Vector of covariate column names
- `random_var`: Column name for random effect grouping (e.g., "patientID")

**Validation checks**:
- Sample ID alignment
- Required variables present
- Factor levels correctly set
- Missing values flagged
- All-zero features identified

### 3. Edge Cases Considered
- **Model convergence failures**: Handled with tryCatch, returns NULL
- **Complete separation**: May occur with sparse PA data
- **Unbalanced designs**: Different timepoints across cohorts
- **No significant results**: Returns empty plot with message
- **Large feature sets**: May need parallelization (future)

---

## Next Steps (TODO)

### Immediate Tasks

1. **Create Helper Functions for Contrasts** (`R/contrasts.R`)
   - `extract_non_reference_contrasts()`: Extract non-reference group trajectory
   - `extract_between_group_contrasts()`: Extract group differences per timepoint
   - Based on existing emmeans code from original script

2. **Create Main Fitting Function** (`R/fit_longitudinal_PA.R`)
   ```r
   fit_longitudinal_PA(
     PA_matrix, metadata, dx_var, timepoint_var,
     covariates, random_var, baseline_level = "t0",
     adjust_method = "BH", return_models = FALSE
   )
   ```
   - Returns list with: `$main_model`, `$non_reference`, `$between_group`, `$combined`

3. **Create Plotting Functions** (`R/plotting.R`)
   - Port `auto_heatmaps_detected()` from original script
   - Port `make_family_heatmap()`
   - Add customization options (colors, themes, filters)

4. **Add Roxygen2 Documentation**
   - Document all exported functions
   - Add examples
   - Generate .Rd files in `man/`

5. **Create Example/Vignette**
   - Simulated example data
   - Step-by-step workflow
   - Interpretation guide

6. **Testing**
   - Unit tests for validation functions
   - Integration tests with example data
   - Edge case tests (convergence failures, etc.)

### Future Enhancements

- **Abundance models**: `fit_longitudinal_abundance()` for count data
- **PA matrix creation**: `create_PA_matrix()` wrapper for CPM-based thresholding
- **Parallel processing**: For large feature sets
- **Model diagnostics**: Residual plots, convergence checks
- **Export functions**: Save results to Excel/CSV with formatting
- **Interactive visualizations**: Plotly integration

---

## File Organization Plan

```
R/
├── utils.R              ✅ Created - validation and helpers
├── contrasts.R          ⏳ Next - emmeans contrast extraction
├── fit_PA.R             ⏳ Main PA model fitting
├── fit_abundance.R      ⏳ Abundance model fitting
├── plotting.R           ⏳ Heatmap functions
└── create_PA_matrix.R   ⏳ CPM-based PA matrix creation

man/                     ⏳ To be generated by roxygen2

tests/testthat/          ⏳ To be created
├── test-validation.R
├── test-contrasts.R
└── test-fitting.R

vignettes/               ⏳ Future
└── tutorial.Rmd
```

---

## Example Usage (Planned)

```r
library(LongitudinalSmarter)

# 1. Fit comprehensive model
results <- fit_longitudinal_PA(
  PA_matrix = total.PA,
  metadata = total.contig.metadata.clean_0.75_0.03,
  dx_var = "Dx.Status",
  timepoint_var = "onset_timeline_combined",
  covariates = c("Country", "Sex", "Age.at.Gluten.Introduction..months.",
                 "HLA.Category", "feeding_first_year", "Delivery.Mode"),
  random_var = "patientID",
  baseline_level = "t0",
  adjust_method = "BH"
)

# Access different result components
results$main_model          # Reference group + interactions
results$non_reference       # Non-reference group trajectory
results$between_group       # Group differences at each timepoint
results$combined            # All combined

# 2. Filter significant results
sig_results <- results$combined %>%
  filter(p.adj < 0.05)

# 3. Generate heatmaps
plots <- plot_longitudinal_heatmaps(
  results = results$combined,
  term_families = c("TP (CONTROL)", "TP (CELIAC)", "Between groups"),
  fdr_threshold = 0.05,
  value_type = "log_odds"
)

# View individual plots
plots$`TP (CONTROL)`
plots$`TP (CELIAC)`
plots$`Between groups`
```

---

## Notes and Decisions

### Why emmeans?
- Allows extraction of non-reference group effects not directly in model output
- Provides clean contrasts with proper SEs and CIs
- Handles complex marginal means calculations

### Why separate non-reference and between-group chunks?
- **Non-reference**: Within-group trajectory (CELIAC vs its own baseline)
- **Between-group**: Cross-group comparison at each timepoint
- Different scientific questions, different FDR correction strategies

### FDR Correction Strategy
- Original script: Group by `term` for main model
- Non-reference contrasts: Group by `term` (e.g., all "TP (CELIAC): t0-6 vs t0" together)
- Between-group: Group by `timepoint` (all contigs at t0-6 together)

### Why Not Just Use DESeq2/edgeR Directly?
- Need mixed-effects for repeated measures (random effect = patientID)
- Need flexible covariate adjustment
- Need specific contrasts not easily obtained from standard RNA-seq tools
- glmmTMB handles sparse PA data better than GLM-based approaches

---

## Timeline

- **2025-01-04**: Project initiated, package structure created
- **Next session**: Complete contrast extraction functions
- **Goal**: v0.1.0 release with core PA functionality

---

## Questions to Resolve

1. Should we store full model objects or just tidied results? (Storage vs flexibility)
2. Parallel processing: Should it be default or opt-in?
3. How to handle cohort-specific analyses (US, Italy, Total)?
4. Should we include example data in the package?
5. CRAN submission plan?

---

## Contact & Contribution

Package author: [To be filled]
Repository: [To be created on GitHub]
Issues: [GitHub Issues]
