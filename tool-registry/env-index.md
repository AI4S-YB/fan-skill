# Environmental Index Construction (环境指数构建)

**Goal:** Transform raw climate/soil variables into interpretable environmental indices
**Best for:** Multi-environment trials with environmental covariate data

## Prerequisites
- R 4.0+, packages: `FactoMineR`, `factoextra`, `stats`, `tidyverse`
- Environmental covariate table: rows = environments, columns = climate/soil variables
- Variables should be numeric and complete (no missing values)

## Basic Usage

### PCA-based Composite Index (≥3 climate vars)

```r
library(FactoMineR)
library(factoextra)
library(tidyverse)

# Load environmental covariate data
env_data <- read.csv("env_covariates.csv", row.names = 1)
# Rows: environments (e.g., "LOC_YEAR")
# Columns: climate/soil variables (e.g., "GDD", "Tmax_mean", "Precip_sum")

# Scale variables (mean=0, sd=1)
env_scaled <- scale(env_data)

# Run PCA
pca_res <- PCA(env_scaled, scale.unit = FALSE, graph = FALSE)

# Inspect variance explained
pca_res$eig  # Eigenvalues and % variance

# Extract PC scores as environmental indices
env_index <- as.data.frame(pca_res$ind$coord[, 1:3])
colnames(env_index) <- c("PC1", "PC2", "PC3")
env_index$Environment <- rownames(env_data)

# PC1 typically = moisture gradient
# PC2 typically = temperature gradient

# Check: PC1+PC2 should explain > 60% variance
var_explained <- sum(pca_res$eig[1:2, 2])
cat(sprintf("PC1+PC2 variance explained: %.1f%%\n", var_explained))

# Variable contributions to each PC
pca_res$var$contrib[, 1:3]

# Save environmental index
write.csv(env_index, "env_index.csv", row.names = FALSE)
```

### Direct Variable Usage (< 3 climate vars)

```r
# When you have only 1-2 climate variables, skip PCA.
# Use raw or scaled variables directly as environmental descriptors.

# Example: only GDD and Precip available
env_index <- env_data %>%
  select(Environment, GDD, Precip_sum) %>%
  mutate(
    GDD_scaled = scale(GDD)[, 1],
    Precip_scaled = scale(Precip_sum)[, 1]
  )

# Use GDD_scaled and Precip_scaled independently
# in reaction norm or regression models.
```

## Plant Relevance

- **Crop growth windows matter**: Do NOT use whole-season averages blindly. Define crop-specific critical windows (e.g., flowering period for cereals, grain-filling for legumes).
- **Variable groups**: Group variables by biological function before PCA:
  - Temperature group (GDD, Tmean, Tmax, Tmin, diurnal range)
  - Moisture group (total precipitation, rainy days, VPD, soil moisture)
  - Radiation group (solar radiation, sunshine hours, PAR)
- **Separate PCAs per group** can be more interpretable than one mega-PCA.
- **Common climate variables** in plant breeding:
  - W1201 (Tmax during vegetative phase), W1202 (Tmin during vegetative)
  - W1203 (precipitation during reproductive phase)
  - GDD (growing degree days, base temperature crop-specific)
  - VPD (vapor pressure deficit) during critical period

## When to Use Each Method

| Scenario | Method | Reason |
|----------|--------|--------|
| ≥3 climate variables | PCA | Dimension reduction, captures collinearity |
| 1-2 climate variables | Direct use | PCA not meaningful |
| Many variables + clear groups | Grouped PCA | Biological interpretability |
| Time-series climate data | Functional PCA | Captures temporal dynamics |

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| PC1 explains > 80% | One variable dominates (e.g., elevation) | Check variable scaling; consider removing extreme variables |
| PC1+PC2 < 50% | Variables are uncorrelated | Use individual variables instead; PCA not helpful |
| Missing values in input | Incomplete climate records | Impute from nearby stations or reanalysis data |
| PCA sample size < variables | More variables than environments | Reduce variables or use regularized PCA |
