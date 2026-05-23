# BLUP / BLUE Estimation

**Goal:** Estimate best linear unbiased predictions (BLUP) or best linear unbiased estimates (BLUE) for genotypes
**Best for:** Ranking genotypes, GS training data preparation, multi-environment trials

## Prerequisites
- R 4.0+, lme4 or sommer
- Phenotype data with genotype identifiers

## Multi-Environment BLUP

```r
library(lme4)
blup_model <- lmer(trait ~ environment + (1|genotype) + (1|environment:genotype), data = pheno)
blups <- ranef(blup_model)$genotype
```

## Single-Environment BLUP with sommer

```r
library(sommer)
ans <- mmer(trait ~ 1, random = ~ vs(genotype), data = pheno)
blups <- ans$U[[1]]$trait
```

## BLUE (Fixed Genotype)

```r
blue_model <- lm(trait ~ genotype + environment, data = pheno)
blues <- coef(blue_model)
```

## Plant-Specific Notes
- BLUP shrinks estimates toward the mean — better for ranking than absolute values
- For GS training data, use BLUPs across environments as the response variable
- Reliability (r²) of BLUPs should be reported alongside estimates

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "Model failed to converge" | Unbalanced design or insufficient data | Simplify random effects structure |
