# Stage/Tissue Specificity: Tau Index

**Goal:** Quantify how specifically a gene is expressed across tissues
or developmental stages.

**Best for:** Tissue atlases, organ comparisons, developmental stage
series with >= 3 conditions.

## Prerequisites

- R 4.0+, base R (no additional packages required)
- Expression matrix: genes x tissues/stages, normalized (TPM, FPKM, or VST)
- At least 3 conditions for meaningful tau calculation

## Basic Usage

```r
# Tau index: 0 = ubiquitous, 1 = strictly tissue-specific
# Formula: tau = sum(1 - (xi / max(x))) / (n - 1)

tau_index <- function(expr_vector) {
  x <- expr_vector
  max_x <- max(x)

  if (max_x == 0) return(NA)  # no expression anywhere

  x_norm <- x / max_x          # normalize to max = 1
  tau <- sum(1 - x_norm) / (length(x) - 1)

  return(tau)
}

# Apply to expression matrix (rows = genes, cols = tissues)
tau_values <- apply(expr_matrix, 1, tau_index)

# Remove NA genes (no expression)
tau_values <- tau_values[!is.na(tau_values)]

# Summary
cat("Number of genes analyzed:", length(tau_values), "\n")
cat("Median tau:", median(tau_values), "\n")
cat("Tissue-specific (tau > 0.8):", sum(tau_values > 0.8), "genes\n")
cat("Tissue-enriched (tau 0.5-0.8):",
    sum(tau_values > 0.5 & tau_values <= 0.8), "genes\n")
cat("Broadly expressed (tau < 0.5):", sum(tau_values < 0.5), "genes\n")

# Histogram
pdf("tau_distribution.pdf", width = 7, height = 5)
hist(tau_values, breaks = 50, col = "steelblue", border = "white",
     xlab = "Tau index", ylab = "Number of genes",
     main = "Expression Specificity Distribution")
abline(v = c(0.5, 0.8), lty = 2, col = c("orange", "red"))
legend("topright", legend = c("tau = 0.5 (enriched)", "tau = 0.8 (specific)"),
       lty = 2, col = c("orange", "red"))
dev.off()

# Identify the tissue with maximum expression for each gene
max_tissue <- apply(expr_matrix, 1, function(row) {
  idx <- which.max(row)
  if (length(idx) == 0) return(NA)
  return(colnames(expr_matrix)[idx])
})

# Tissue-specific genes + which tissue
ts_genes <- data.frame(
  gene = names(tau_values[tau_values > 0.8]),
  tau  = tau_values[tau_values > 0.8],
  max_tissue = max_tissue[tau_values > 0.8],
  stringsAsFactors = FALSE
)
head(ts_genes[order(-ts_genes$tau), ], 20)

# Tissue-specificity enrichment per tissue
tissue_counts <- table(max_tissue)
tissue_ts_genes <- sapply(names(tissue_counts), function(tis) {
  sum(tau_values[max_tissue == tis] > 0.8, na.rm = TRUE)
})
print(data.frame(
  tissue = names(tissue_counts),
  total_genes = as.integer(tissue_counts),
  tissue_specific = as.integer(tissue_ts_genes),
  pct_specific = round(100 * tissue_ts_genes / tissue_counts, 1)
))
```

## Tissue-Specificity Score (TSI) -- Alternative

```r
# TSI = max(xi) / sum(xi), alternative to tau
# Range: 1/n (ubiquitous) to 1 (specific)
tsi <- function(expr_vector) {
  x <- expr_vector
  if (sum(x) == 0) return(NA)
  return(max(x) / sum(x))
}

tsi_values <- apply(expr_matrix, 1, tsi)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| tau threshold (specific) | 0.8 | Genes with tau > this are "tissue-specific" |
| tau threshold (enriched) | 0.5 | Genes with tau > 0.5 are "tissue-enriched" |
| minimum expression | > 0 | Genes with zero max expression are excluded |

## Plant-Specific Notes

- **Tissue atlases** (e.g., rice eFP, Arabidopsis TraVA, maize Atlas):
  tau is the standard metric. TSI is more interpretable when you want
  "what fraction of total expression is in the top tissue."
- **Developmental series** (e.g., leaf 1-8, embryo stages): tau
  measures temporal specificity. Genes with high tau in a developmental
  series are "stage-specific" (often TFs).
- **Caveat -- tissue set dependency**: tau depends on which tissues
  are included. A "leaf-specific" gene in a leaf-root-flower comparison
  may appear broadly expressed when seed and meristem are added.
  Always report the tissue panel used.
- **Polyploid crops**: calculate tau separately per homeolog and
  compare. Homeolog expression bias (one copy more tissue-specific
  than the other) is common in polyploid plants.
- **Photosynthetic genes**: typically low tau (expressed in all
  green tissues). High-tau photosynthetic genes are suspicious --
  check if they are actually C4-specific or bundle-sheath-specific.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| All tau values near 0 | Too few tissues included | Minimum 3 tissues; ideally 5+ |
| All tau values near 1 | Many genes expressed in only one tissue | Check data quality; confirm normalization |
| `NA` for many genes | Zero expression across all tissues | Filter genes with max expression > threshold |
| Tau > 1 | Calculation error | Check denominator: n_tissues - 1, not n_tissues |
| Tissue-specific genes concentrated in one tissue | Normalization bias for one tissue | Check that per-tissue median expression is similar |
