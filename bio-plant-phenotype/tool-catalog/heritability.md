# Heritability Estimation

**Goal:** Estimate broad-sense (H²) or narrow-sense (h²) heritability from phenotype data
**Best for:** All quantitative traits — essential before GWAS or GS

## Prerequisites
- R 4.0+, lme4, sommer
- Phenotype data with genotype (or line/variety) identifiers and replication structure

## Mixed Model H² (Multi-Environment)

```r
library(lme4)
# Broad-sense heritability with genotype as random effect
model <- lmer(trait ~ (1|genotype) + (1|environment) + (1|environment:rep), data = pheno)
vc <- as.data.frame(VarCorr(model))
H2 <- vc$vcov[1] / sum(vc$vcov)
```

## ANOVA-based H² (Single Environment)

```r
# One-way ANOVA
aov_model <- aov(trait ~ genotype, data = pheno)
ms <- summary(aov_model)[[1]][, "Mean Sq"]
H2 <- (ms[1] - ms[2]) / (ms[1] + (n_reps - 1) * ms[2])
```

## Plant-Specific Notes
- Self-pollinated crops: genetic variance among inbred lines is the primary component
- Multi-environment trials provide more reliable H² estimates
- Single-environment H² is environment-specific and not generalizable

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "H² > 1" | Negative variance component or unbalanced design | Check design balance |
| "H² ≈ 0" | Low genetic variance or high error | May be real — check CV |
