# edgeR

**Goal:** Differential expression using empirical Bayes quasi-likelihood
**Best for:** Low-replicate designs (2 replicates), or when TMM normalization is preferred

## Prerequisites
- R 4.0+, edgeR (Bioconductor)

## Basic Usage

```r
library(edgeR)

# Create DGEList
y <- DGEList(counts = counts, group = group)

# Filter low-expression genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize (TMM)
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateDisp(y)

# Test (quasi-likelihood F-test)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, contrast = contrast)
topTags(qlf, n = 100)
```

## When to Choose edgeR over DESeq2

- **2 replicates**: edgeR's QL framework may perform slightly better
- **No replicates**: edgeR can estimate dispersion from the data (with strong caveats)
- **TMM normalization preferred**: Some plant datasets with strong compositional bias benefit from TMM

## Plant-Specific Notes

- edgeR's TMM normalization can be beneficial for plant samples with very different RNA composition (e.g., seeds vs leaves)
- Works with the same count matrix format as DESeq2

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "no residual df" | No replicates | Use `estimateGLMCommonDisp()` with caution |
| "contrast vector wrong length" | Design matrix mismatch | Check `colnames(design)` |
