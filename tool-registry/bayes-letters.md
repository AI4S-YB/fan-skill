# BayesA / BayesB / BayesC (via BGLR)

**Goal:** Bayesian genomic prediction with variable selection
**Best for:** High-density markers (≥50K), moderate-large populations

## Prerequisites
- R 4.0+, BGLR package
- Numeric genotype matrix
- Continuous phenotype

## Basic Usage (BayesB)

```r
library(BGLR)

# BayesB: mixture prior — some markers have effect, some don't
ETA <- list(list(X = geno_matrix, model = "BayesB"))

fm <- BGLR(
  y = pheno,
  ETA = ETA,
  nIter = 20000,
  burnIn = 5000,
  thin = 10
)

gebv <- fm$ETA[[1]]$b  # marker effects
gv <- fm$yHat           # predicted breeding values
```

## Bayes Model Comparison

| Model | Prior | Best For |
|-------|-------|----------|
| BayesA | Scaled-t (all markers have effect) | Many small-effect QTL |
| BayesB | Mixture (some markers = 0) | Few large-effect QTL |
| BayesC | Mixture with common variance | Intermediate |

## Plant-Specific Notes

- BayesB often outperforms GBLUP by 2-5% when marker density >50K
- For self-pollinated crops (most major QTL fixed), the gain from variable selection is modest
- MCMC convergence should be checked (trace plots) for large plant genomes
