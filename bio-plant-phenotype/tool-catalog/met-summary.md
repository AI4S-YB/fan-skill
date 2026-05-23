# Multi-Environment Trial Summary

**Goal:** Summarize genotype performance across environments, assess stability
**Best for:** Multi-environment trials (≥2 environments)

## AMMI (Additive Main Effects and Multiplicative Interaction)

```r
library(agricolae)
ammi <- AMMI(env, genotype, rep, trait, data = pheno)
plot(ammi)
```

## GGE Biplot

```r
library(GGEBiplots)
gge <- GGEModel(pheno_matrix)
plot(gge, type = 1)  # Which-won-where pattern
plot(gge, type = 2)  # Mean vs stability
```

## Plant-Specific Notes
- AMMI/GGE help identify broadly adapted vs specifically adapted genotypes
- "Which-won-where" biplot shows which genotype performs best in each environment
- For breeding: select genotypes stable across target environments

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "Cannot compute AMMI" | Missing values in phenotype matrix | Impute or remove genotypes with missing data |
