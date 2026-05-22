# limma-voom

**Goal:** Differential expression using linear models with precision weights
**Best for:** Large datasets (≥30 samples), continuous covariates, complex designs

## Prerequisites
- R 4.0+, limma, edgeR

## Basic Usage

```r
library(limma)
library(edgeR)

# Create DGEList and normalize
y <- DGEList(counts = counts)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# Voom transformation
design <- model.matrix(~ 0 + group)
v <- voom(y, design)

# Linear model + empirical Bayes
fit <- lmFit(v, design)
contrast <- makeContrasts(groupTrt - groupCtl, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
topTable(fit2, n = 100)
```

## When to Choose limma-voom

- **>30 samples**: voom is faster than DESeq2 for large datasets
- **Complex designs**: continuous covariates, time series, factorial designs
- **>40,000 genes**: memory efficiency for large genomes

## Plant-Specific Notes

- Good for large plant genomes where the gene count exceeds 40K
- Handles complex breeding designs (diallel, factorial mating) well
