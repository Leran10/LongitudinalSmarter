# Installation and Testing Guide

## ✅ Package Status

The LongitudinalSmarter package is now **ready for GitHub release**!

All critical issues have been fixed:
- ✓ Dependencies added to DESCRIPTION
- ✓ Documentation generated with roxygen2
- ✓ NAMESPACE properly configured
- ✓ Package loads successfully
- ✓ Basic functionality tested and working

---

## Installation from GitHub

Once you push to GitHub, users can install with:

```r
# Install devtools if needed
install.packages("devtools")

# Install LongitudinalSmarter
devtools::install_github("yourusername/LongitudinalSmarter")
```

---

## Local Installation (for testing)

To install and test locally before pushing to GitHub:

```r
# Install from local directory
devtools::install("/Users/leranwang/Documents/github_projects/LongitudinalSmarter")

# Or load without installing (for development)
devtools::load_all("/Users/leranwang/Documents/github_projects/LongitudinalSmarter")
```

---

## Quick Test

After installation, test the package:

```r
library(LongitudinalSmarter)

# Create minimal test data
set.seed(123)
PA_test <- matrix(
  sample(0:1, 200, replace = TRUE),
  nrow = 10, ncol = 20,
  dimnames = list(paste0("feature_", 1:10), paste0("sample_", 1:20))
)

metadata_test <- data.frame(
  Dx.Status = factor(rep(c("CONTROL", "CELIAC"), each = 10)),
  timepoint = factor(rep(c("t0", "t0-6", "t0-12", "t0-18"), 5)),
  patientID = paste0("P", rep(1:10, each = 2)),
  row.names = paste0("sample_", 1:20)
)

# Run analysis
results <- fit_longitudinal_PA(
  PA_matrix = PA_test,
  metadata = metadata_test,
  dx_var = "Dx.Status",
  timepoint_var = "timepoint",
  random_var = "patientID"
)

# Check results
str(results, max.level = 1)
```

---

## Pushing to GitHub

### Step 1: Initialize Git Repository

```bash
cd /Users/leranwang/Documents/github_projects/LongitudinalSmarter

# Initialize git (if not already done)
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: LongitudinalSmarter R package v0.1.0

- Comprehensive longitudinal mixed-effects models for microbiome data
- Automatic extraction of reference and non-reference group trajectories
- Between-group contrasts at each timepoint
- Publication-ready heatmap visualization
- Full documentation and examples included"
```

### Step 2: Create GitHub Repository

1. Go to https://github.com/new
2. Repository name: `LongitudinalSmarter`
3. Description: `R package for comprehensive longitudinal mixed-effects models in microbiome studies`
4. Make it **Public** (so others can install it)
5. **Do NOT** initialize with README (you already have one)
6. Click "Create repository"

### Step 3: Push to GitHub

```bash
# Add the remote
git remote add origin https://github.com/YOURUSERNAME/LongitudinalSmarter.git

# Push to GitHub
git branch -M main
git push -u origin main
```

### Step 4: Update README

After pushing, update the README.md installation section with your actual GitHub username:

```r
devtools::install_github("YOURUSERNAME/LongitudinalSmarter")
```

---

## Files Included

```
LongitudinalSmarter/
├── DESCRIPTION              # Package metadata ✓
├── NAMESPACE               # Exported functions ✓
├── LICENSE                 # MIT license ✓
├── README.md               # User documentation ✓
├── INSTALL.md              # This file ✓
├── claude.md               # Development notes ✓
├── .gitignore              # Git ignore rules ✓
├── R/                      # Source code ✓
│   ├── utils.R
│   ├── contrasts.R
│   ├── fit_PA.R
│   └── plotting.R
├── man/                    # Documentation (10 .Rd files) ✓
├── examples/
│   └── example_usage.R     # Complete workflow example ✓
└── tests/                  # Future unit tests
```

---

## What Users Will See

When users install from GitHub:

1. Install command:
   ```r
   devtools::install_github("YOURUSERNAME/LongitudinalSmarter")
   ```

2. Load and use:
   ```r
   library(LongitudinalSmarter)
   ?fit_longitudinal_PA  # View documentation
   ```

3. Access help:
   ```r
   help(package = "LongitudinalSmarter")
   ```

---

## Next Steps (Optional)

### Recommended before public release:
1. ✓ ~~Fix critical errors~~ (DONE)
2. Test with real data from your celiac study
3. Add unit tests in `tests/testthat/`
4. Create a vignette tutorial
5. Add badges to README (build status, CRAN, etc.)

### Future enhancements:
- Add `fit_longitudinal_abundance()` for count data
- Add `create_PA_matrix()` for CPM-based PA creation
- Parallel processing for large datasets
- More visualization options

---

## Troubleshooting

### "Package not found" error
Make sure the repository is public on GitHub.

### Installation fails
Check that all dependencies are available:
```r
install.packages(c("dplyr", "tidyr", "purrr", "tibble", "stringr",
                   "ggplot2", "glmmTMB", "broom.mixed", "emmeans", "forcats"))
```

### Function not found
Make sure you loaded the package:
```r
library(LongitudinalSmarter)
```

---

## Support

Once on GitHub:
- **Issues**: Open an issue at https://github.com/YOURUSERNAME/LongitudinalSmarter/issues
- **Pull requests**: Contributions welcome!

---

## Citation

```
Wang L (2025). LongitudinalSmarter: Comprehensive Longitudinal
Mixed-Effects Models for Microbiome Data. R package version 0.1.0.
https://github.com/YOURUSERNAME/LongitudinalSmarter
```
