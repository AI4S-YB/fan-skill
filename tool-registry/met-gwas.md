# MET-GWAS (Multi-Environment Trial GWAS)

**Goal:** Joint GWAS across multiple environments — separating main genetic effects from G×E
**Best for:** Multi-environment phenotype data (the norm in plant breeding)

## Approach

Two main strategies:

### Strategy 1: Per-environment + meta-analysis
1. Run GWAS separately for each environment
2. Compare significant SNPs across environments
3. Stable QTL = significant in ≥2 environments
4. Environment-specific QTL = significant in only 1

### Strategy 2: Joint model (BLUP)
Use BLUPs to summarize across environments, then run single GWAS.

```r
# Using GAPIT with BLUP across environments
# Step 1: Calculate BLUP for each genotype across environments
library(lme4)
blup_model <- lmer(trait ~ (1|genotype) + (1|environment), data = pheno)
blups <- ranef(blup_model)$genotype

# Step 2: GWAS on BLUPs
myGAPIT <- GAPIT(G = myGD, Y = blups, model = "CMLM")
```

## Key Decisions

| Scenario | Recommendation |
|----------|---------------|
| 2-3 environments, same management | Per-environment + compare |
| >3 environments, varying conditions | BLUP → GWAS on BLUPs |
| Strong G×E suspected | Per-environment + look for env-specific peaks |
| G×E not a concern | Pool or BLUP |

## Plant-Specific Notes

- Multi-environment trials are standard in plant breeding — always check if data has env info
- Environment can include: location, year, treatment, management practice
- env_count is detected by `inspect_data.sh`
