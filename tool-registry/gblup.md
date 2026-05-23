# GBLUP (Genomic Best Linear Unbiased Prediction)

**Goal:** Genomic prediction using genomic relationship matrix (G matrix)
**Best for:** Most breeding scenarios — the default GS method

## Prerequisites
- R 4.0+, sommer or BGLR package
- Genotype matrix or precomputed G matrix
- Phenotype data (single or multi-environment)

## Single-Trait GBLUP

```r
library(sommer)

# Build G matrix
G <- A.mat(t(geno_matrix) - 1)

# Mixed model
ans <- mmer(
  trait ~ 1,
  random = ~ vs(genotype, Gu = G),
  data = pheno_data
)

gebv <- ans$U[[1]]$trait
```

## Multi-Environment with G×E

```r
# GBLUP with genotype × environment interaction
ans <- mmer(
  trait ~ env,
  random = ~ vs(genotype, Gu = G) + vs(env:genotype),
  data = pheno_data
)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| G matrix | VanRaden method 1 | Standard in plant breeding |
| Multi-trait | `mmer(fixed = cbind(trait1, trait2) ~ 1, ...)` | Uses genetic correlations |

## Plant-Specific Notes

- GBLUP with G×E is the workhorse of modern plant breeding programs
- For multi-environment data, always try the G×E model — it often outperforms single-environment GS
- Balance speed vs accuracy: GBLUP is O(n³) with sample size; for >5000 samples, use AI-REML approximations
