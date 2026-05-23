# Spatial Correction for Field Trials

**Goal:** Correct for spatial trends in field trials using row-column models
**Best for:** Field trials with known row/column positions, CV > 15%

## Prerequisites
- R 4.0+, sommer or SpATS
- Row and column coordinates for each plot

## AR1×AR1 Spatial Model (sommer)

```r
library(sommer)
ans <- mmer(trait ~ genotype,
  random = ~ vs(row) + vs(col) + vs(units),
  rcov = ~ vs(units),
  data = pheno)
```

## SpATS (Spatial Analysis of Field Trials with Splines)

```r
library(SpATS)
ans <- SpATS(response = "trait", genotype = "genotype",
  spatial = ~ PSANOVA(row, col, nseg = c(10, 10)),
  data = pheno)
```

## Plant-Specific Notes
- Spatial variation is common in large breeding trials — always check CV
- Row-column designs (lattice, alpha) already partially control for spatial effects
- Spatial correction can substantially increase heritability estimates

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "Model did not converge" | Too few rows/columns or extreme imbalance | Use spline-based SpATS instead |
