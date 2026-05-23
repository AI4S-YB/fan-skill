# Multi-Tissue eQTL Analysis

**Goal:** Classify eQTLs as tissue-specific or tissue-shared, quantify sharing
patterns, and identify tissue-specific regulatory mechanisms.

**Best for:** Studies with expression data from >= 2 tissues, especially plant
developmental atlases or stress-response panels.

## Prerequisites

- R 4.0+, mashr, metaCCA, or Meta-Tissue packages
- eQTL summary statistics per tissue (beta + SE for all SNP-gene pairs)
- Consistent sample overlap information between tissues
- Tissue metadata (names, categories: vegetative/reproductive/stress)

## Basic Usage: mashR for Multi-Tissue eQTL

```r
library(mashr)
library(ashr)

# --- Input: matrix of effects and standard errors per tissue ---
# Each row = one eQTL (SNP-gene pair in at least one tissue)
# Each column = one tissue pair (beta + SE alternating)

# Example: 3 tissues -> 6 columns (beta_T1, SE_T1, beta_T2, SE_T2, beta_T3, SE_T3)
# Load from per-tissue MatrixEQTL outputs

tissues <- c("leaf", "root", "seed")
n_tissues <- length(tissues)

# Build mash data object
beta_matrix <- as.matrix(eqtl_combined[, grep("^beta_", colnames(eqtl_combined))])
se_matrix   <- as.matrix(eqtl_combined[, grep("^SE_", colnames(eqtl_combined))])

data <- mash_set_data(beta_matrix, se_matrix)

# --- Estimate data-driven covariance matrices ---
# 1. Identify strong eQTLs (most significant per condition)
strong_subset <- get_significant_results(data, thresh = 0.05)

# 2. Compute PCA-based covariance patterns
data_pca <- mash_set_data(
  beta_matrix[strong_subset, ],
  se_matrix[strong_subset, ]
)
U_pca <- cov_pca(data_pca, n_tissues)

# 3. Extreme deconvolution for additional patterns
U_ed <- cov_ed(data, U_pca, strong_subset)

# --- Fit mash model ---
U_canonical <- cov_canonical(data)

m <- mash(data, c(U_canonical, U_ed))

# --- Tissue-specific vs shared classification ---
# Extract posterior summaries
pm <- get_pm(m)   # posterior means
psd <- get_psd(m) # posterior standard deviations
lfsr <- get_lfsr(m) # local false sign rate

# For each eQTL, compute sharing statistic:
# 1. Significant in how many tissues? (lfsr < 0.05)
sharing <- apply(lfsr, 1, function(x) sum(x < 0.05))

eqtl_classification <- data.frame(
  eqtl_id = eqtl_combined$eqtl_id,
  n_tissues_significant = sharing,
  type = ifelse(sharing == n_tissues, "shared",
         ifelse(sharing == 1, "tissue_specific", "partially_shared"))
)

table(eqtl_classification$type)
```

## Alternative: Meta-Analysis Approach

For simpler cases where mashR is overkill:

```r
library(metafor)

# Fisher's method to combine p-values across tissues
combine_pvalues <- function(pvalues) {
  # Fisher's method
  chi_sq <- -2 * sum(log(pvalues))
  df <- 2 * length(pvalues)
  pchisq(chi_sq, df, lower.tail = FALSE)
}

# Per eQTL: combine tissue-level p-values
eqtl_combined$meta_pvalue <- apply(
  eqtl_combined[, grep("^pvalue_", colnames(eqtl_combined))],
  1, combine_pvalues
)

# FDR correction
eqtl_combined$meta_fdr <- p.adjust(eqtl_combined$meta_pvalue, method = "BH")

# Tissue-specificity score (I-squared heterogeneity)
calc_tissue_specificity <- function(betas, ses) {
  if (length(betas) <= 1) return(0)
  # Random-effects meta-analysis
  res <- rma(yi = betas, sei = ses, method = "REML")
  return(res$I2)  # I^2 = 0% (shared) to 100% (tissue-specific)
}

eqtl_combined$I2 <- apply(eqtl_combined, 1, function(row) {
  betas <- as.numeric(row[grep("^beta_", colnames(eqtl_combined))])
  ses   <- as.numeric(row[grep("^SE_", colnames(eqtl_combined))])
  calc_tissue_specificity(betas, ses)
})
```

## Tissue-Sharing Visualization

```r
library(UpSetR)
library(pheatmap)

# UpSet plot: which tissue combinations share eQTLs?
tissue_matrix <- ifelse(lfsr < 0.05, 1, 0)
colnames(tissue_matrix) <- tissues
upset(as.data.frame(tissue_matrix), nsets = n_tissues,
      order.by = "freq", main.bar.color = "#2166AC")

# Heatmap: correlation of eQTL effect sizes across tissues
cor_matrix <- cor(get_pm(m), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE,
         main = "eQTL Effect Size Correlation Across Tissues",
         color = colorRampPalette(c("#B2182B", "white", "#2166AC"))(100))
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| n_tissues minimum | 2 | Even 2 tissues enables shared/specific classification |
| mashR: strong subset threshold | p < 0.05 in >= 1 tissue | Select informative eQTLs for covariance learning |
| lfsr threshold | < 0.05 | Local false sign rate for significance |
| I^2 threshold (tissue-specific) | >= 75% | High heterogeneity = tissue-specific |
| I^2 threshold (shared) | < 25% | Low heterogeneity = tissue-shared |
| Min samples per tissue | >= 30 | Power drops sharply below 30 |

## Plant-Specific Notes

- **Polyploid multi-tissue**: Run per subgenome per tissue. Homeolog expression
  patterns differ by tissue — a shared eQTL in subgenome A may have a tissue-specific
  homeolog counterpart in subgenome B.
- **Developmental tissue series**: Plant tissues are often developmental stages
  (e.g., seedling -> tillering -> heading -> grain fill). These are ordered —
  consider trend analysis in addition to pairwise sharing.
- **Stress vs control pairs**: The most informative plant multi-tissue design is
  control vs stress for the same tissue. eQTLs appearing only under stress are
  candidates for stress-responsive regulatory variation.
- **Diurnal sampling**: If tissues were collected at different times of day, diurnal
  expression patterns may inflate tissue-specific eQTL calls. Account for sampling
  time as a covariate.
- **Tissue atlas depth**: Many plant atlases are broad (>20 tissues) but shallow
  (n=2-3 per tissue). Multi-tissue eQTL analysis with such data is exploratory only.
- **mashR convergence**: mashR requires substantial compute for full eQTL datasets.
  For > 100K eQTLs across > 5 tissues, use the random subset approach for covariance
  estimation, then apply to the full set.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| mashR: "system is computationally singular" | Too few strong signals or collinear tissues | Check tissue correlation; increase strong subset threshold |
| I^2 = NA for single-tissue eQTLs | Only detected in one tissue | Assign as "tissue-specific" by definition |
| All eQTLs classified as tissue-specific | Low sample overlap between tissues | Verify that the same individuals were used across tissues |
| mashR memory error | Full dataset too large | Use random subset (10K eQTLs) for covariance learning |
| Zero shared eQTLs across unrelated tissues | Biologically expected (e.g., leaf vs root) | Focus on within-organ comparisons |
