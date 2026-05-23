# GAPIT FarmCPU (Fixed and Random Model Circulating Probability Unification)

**Goal:** GWAS with iterative fixed + random model — high power, low false positive
**Best for:** Medium-to-large populations with population stratification (λ > 1.05)

## Prerequisites
- R 4.0+, GAPIT 3.0+
- Marker density ≥ 1000 SNPs (requires enough markers for polygenic background)

## Basic Usage

```r
library(GAPIT3)

myGAPIT <- GAPIT(
  G = myGD,
  GD = myGD,
  GM = myGM,
  Y = myY,
  model = "FarmCPU",
  PCA.total = 3,
  Multiple_analysis = TRUE   # Run multiple iterations
)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| PCA.total | 3-5 | Controls population structure |
| Multiple_analysis | TRUE | Multiple iterations improve stability |

## When FarmCPU Shines

- Population stratification is present but you want high statistical power
- Marker density is adequate (>1000 SNPs) — FarmCPU needs markers to estimate polygenic effects
- You want fewer false positives than GLM but more power than strict CMLM

## When to Avoid

- **Low marker density (<1000 SNPs)**: FarmCPU can't properly model the polygenic background
- **Very small samples (<100)**: Model fitting unstable

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Model did not converge" | Insufficient markers or samples | Switch to CMLM |
| "Too few markers after filtering" | QC removed too many markers | Relax QC thresholds |
