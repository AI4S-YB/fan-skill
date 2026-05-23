# Differential Protein Abundance Analysis

**Goal:** Identify significantly changed proteins between experimental conditions
**Best for:** Plant proteomics data with >= 3 biological replicates per condition

## Prerequisites

- Protein quantification matrix (from MaxQuant, DIA-NN, or Spectronaut)
- Sample metadata table (conditions, batches)
- R 4.0+ with limma, Biobase

## Basic Usage (limma)

```r
library(limma)
library(tidyverse)

# Step 1: Load and prepare data
lfq <- read.delim("proteinGroups.txt", stringsAsFactors = FALSE)

# Filter: remove contaminants, reverse hits, "only identified by site"
lfq_clean <- lfq %>%
  filter(Reverse != "+",
         Potential.contaminant != "+",
         Only.identified.by.site != "+")

# Extract LFQ intensities
lfq_cols <- grep("^LFQ.intensity", colnames(lfq_clean), value = TRUE)
lfq_mat <- lfq_clean[, lfq_cols] %>%
  as.matrix() %>%
  log2()
lfq_mat[is.infinite(lfq_mat)] <- NA

# Step 2: Filter for valid values
# Require >= 2 valid values in at least one condition
valid_count <- rowSums(!is.na(lfq_mat))
lfq_filt <- lfq_mat[valid_count >= 3, ]

# Step 3: Normalize (median normalization)
lfq_norm <- normalizeBetweenArrays(lfq_filt, method = "quantile")

# Step 4: Missing value imputation (KNN)
library(impute)
lfq_imp <- impute.knn(lfq_norm, k = 10)$data

# Step 5: Differential analysis
# Define groups
groups <- factor(c("Control", "Control", "Control",
                   "Treatment", "Treatment", "Treatment"))

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(lfq_imp, design)

# Contrast
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Results
results <- topTable(fit2, coef = 1, number = Inf,
                    adjust.method = "BH", sort.by = "P")
results$significant <- results$adj.P.Val < 0.05 & abs(results$logFC) > 1

# Summary
table(results$significant)
```

## Missing Value Strategy

```r
# Classify missing values
na_pattern <- rowSums(is.na(lfq_filt)) / ncol(lfq_filt)

# MAR: Left-censored imputation (MinProb)
# MNAR: Remove from analysis or impute with low values

# MinProb imputation for MNAR
impute_minprob <- function(mat, q = 0.01) {
  for (i in seq_len(nrow(mat))) {
    row_vals <- mat[i, ][!is.na(mat[i, ])]
    if (length(row_vals) == 0) next
    min_val <- quantile(row_vals, q, na.rm = TRUE)
    shrink <- min_val - 0.3 * sd(row_vals, na.rm = TRUE)
    mat[i, ][is.na(mat[i, ])] <- rnorm(sum(is.na(mat[i, ])),
                                        mean = shrink, sd = 0.3 * sd(row_vals, na.rm = TRUE))
  }
  return(mat)
}
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| P-value threshold (adj) | 0.05 | Standard FDR |
| |log2FC| threshold | 1 (2-fold) | Standard for abundant proteins |
| Missing value filter | >= 3 valid values in at least 1 condition | Balance sensitivity and specificity |
| Imputation k (KNN) | 10 | Trade-off between accuracy and stability |
| Normalization | quantile or median | Depending on data distribution |

## Plant-Specific Notes

- For phosphoproteomics: |log2FC| >= 1 is too strict; use >= 0.58 (1.5-fold)
- Proteins with only one sample having valid values often represent noise, not biology
- Check if differentially abundant proteins align with RuBisCO depletion efficiency
- Compare DE proteins to RNA-seq DEGs: low correlation (R^2 ~0.3-0.5) is expected and normal
- Plant secondary metabolites may cause run-specific missing values — always visualize by run

## QC Visualization

```r
# PCA plot
pca <- prcomp(t(lfq_imp), scale. = TRUE)
plot(pca$x[, 1:2], col = groups, pch = 19, cex = 2)
legend("topright", legend = levels(groups), col = 1:2, pch = 19)

# Coefficient of variation by condition
cv_per_condition <- function(mat, groups) {
  cvs <- sapply(levels(groups), function(g) {
    cols <- which(groups == g)
    row_sds <- apply(mat[, cols], 1, sd, na.rm = TRUE)
    row_means <- rowMeans(mat[, cols], na.rm = TRUE)
    row_sds / row_means * 100
  })
  return(cvs)
}

# Check CV distribution
cv_res <- cv_per_condition(lfq_imp, groups)
boxplot(cv_res, ylim = c(0, 100),
        ylab = "Coefficient of Variation (%)")
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No significant proteins" | Too strict threshold, or genuine lack of change | Relax |log2FC| to 0.58, check PCA |
| "design matrix not full rank" | Confounded design | Simplify model, check metadata |
| KNN imputation fails | Too many missing values | Increase valid value filter |
| Unrealistic log2FC | Unnormalized or non-logged data | Check if input is log2 transformed |
