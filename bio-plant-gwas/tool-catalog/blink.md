# BLINK (Bayesian-information and Linkage-disequilibrium Iteratively Nested Keyway)

**Goal:** High-power GWAS with LD-based binning
**Best for:** Large populations (≥500) with high-density markers (≥50K) and low stratification

## Prerequisites
- R 4.0+, GAPIT 3.0+
- High-density markers

## Basic Usage

```r
library(GAPIT3)

myGAPIT <- GAPIT(
  G = myGD,
  GD = myGD,
  GM = myGM,
  Y = myY,
  model = "BLINK",
  PCA.total = 3
)
```

## When BLINK is Best

- Large sample size (≥500)
- High marker density (≥50K SNPs)
- Low population stratification (λ < 1.05)
- You want maximum statistical power

## When to Avoid

- Low marker density: BLINK's LD-based binning requires many markers
- Strong population stratification: BLINK's default PCA correction may be insufficient — use FarmCPU instead
- Small samples: parameter estimation unstable

## Plant Relevance

- Good for: modern breeding populations (large, high-density genotyping), maize NAM panels
- Not ideal for: small biparental populations, low-density legacy datasets
