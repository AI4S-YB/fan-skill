# rrBLUP (Ridge Regression BLUP)

**Goal:** Genomic prediction using ridge regression — simplest GS model
**Best for:** Small populations (<200), all marker densities

## Prerequisites
- R 4.0+, rrBLUP package
- Numeric genotype matrix (samples × markers, coded 0/1/2 or -1/0/1)
- Phenotype vector

## Basic Usage

```r
library(rrBLUP)

# Build genomic relationship matrix
G <- A.mat(t(geno_matrix) - 1)

# Mixed model: trait = mean + genotype (random, with G matrix)
ans <- mixed.solve(y = pheno, K = G)

# GEBV = genomic estimated breeding values
gebv <- ans$u
accuracy <- cor(gebv, pheno)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| G matrix build | `A.mat()` with -1 centering | Standard for rrBLUP |
| Alternative | `kin.blup()` | Direct BLUP without explicit G |

## When to Use

- Small training populations (<200) — rrBLUP is the most robust choice
- Low marker density — doesn't rely on variable selection
- Quick baseline accuracy estimate before trying more complex methods

## Plant-Specific Notes

- For inbred crops, use additive relationship matrix (A.mat) — dominance rarely matters
- rrBLUP and GBLUP are mathematically equivalent when using the same G matrix
- Good for initial GS feasibility assessment with limited data

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "system is computationally singular" | Highly correlated markers | LD-prune first or use GBLUP |
| "non-conformable arguments" | Genotype/phenotype sample mismatch | Check sample IDs |
