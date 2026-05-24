# BayesA / BayesB / BayesC — Bayesian Genomic Prediction

**Goal:** Bayesian genomic prediction with variable selection — markers with small effects are shrunk to zero, capturing major QTL
**Approach:** MCMC sampling → posterior distribution of marker effects → GEBV = sum of marker effects × genotypes
**Best for:** High-density markers (>=50K), moderate-to-large populations, traits with known major QTL

## Prerequisites
- R 4.0+, BGLR package
- Numeric genotype matrix (samples x markers, centered)
- Continuous phenotype vector
- Patience: MCMC chains are computationally intensive for plant genomes (>100K markers)

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(BGLR)

# ── Step 1: Prepare data ──
# geno: samples × markers, center/scale recommended
# pheno: numeric vector, no missing values in training set
geno <- ${GENO_MATRIX}
pheno <- ${PHENO_VECTOR}

# Handle missing phenotypes (prediction targets)
y_na <- which(is.na(pheno))
y_obs <- which(!is.na(pheno))

# ── Step 2: Define predictor list (ETA) ──
# BayesB: mixture prior — most markers have zero effect, some have non-zero
ETA <- list(
  list(
    X = geno,
    model = "${BAYES_MODEL}",     # "BayesB", "BayesC", or "BayesA"
    saveEffects = TRUE            # keep marker effects for post-hoc analysis
  )
)

# ── Step 3: Run MCMC ──
fm <- BGLR(
  y = pheno,
  ETA = ETA,
  nIter = ${N_ITER},              # total MCMC iterations, e.g. 20000
  burnIn = ${BURN_IN},            # discarded iterations, e.g. 5000
  thin = ${THIN},                 # thinning interval, e.g. 10
  saveAt = "${OUTPUT_PREFIX}",     # optional: save posterior samples to disk
  verbose = TRUE                  # monitor progress (recommended)
)

# ── Step 4: Extract results ──
# Predicted breeding values (genomic values)
gebv <- fm$yHat

# Marker effects (posterior means)
marker_effects <- fm$ETA[[1]]$b

# Posterior standard deviations of marker effects
marker_sd <- fm$ETA[[1]]$SD.b

# Residual variance
residual_var <- fm$varE

# Genetic variance (from marker effects)
genetic_var <- fm$ETA[[1]]$varB

# ── Step 5: Identify top markers ──
# BayesB: proportion of iterations where marker had non-zero effect
# (stored in fm$ETA[[1]]$d when available, depending on model)
# Alternative: use abs(marker_effects) / marker_sd as signal-to-noise ratio
top_markers <- order(abs(marker_effects), decreasing = TRUE)[1:${N_TOP}]

# ── Step 6: Accuracy assessment ──
if (length(y_na) > 0) {
  # Prediction accuracy (if validation set with NA phenotypes was included)
  accuracy <- cor(gebv[y_na], ${TRUE_VALUES}[y_na])
} else {
  # Cross-validation must be run in separate BGLR calls (see CV strategy below)
  cat("Run k-fold CV with separate BGLR calls for accuracy estimate.\n")
}
```

---

## Bayesian Model Selection

| Model | Prior on marker effects | pi parameter | Best for |
|-------|------|:---:|----------|
| **BayesA** | Scaled-t (all markers have effect, variance differs) | N/A — all markers included | Polygenic traits, many small-effect QTL |
| **BayesB** | Mixture: point-mass at zero + scaled-t slab | \(\pi \approx 0.95\text{--}0.99\) | Few large-effect QTL (oligogenic traits) |
| **BayesC** | Mixture: point-mass at zero + normal slab (common variance) | \(\pi \approx 0.95\text{--}0.99\) | Intermediate; more stable than BayesB |

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `nIter` | 20000 | <10K markers → 10000; >100K markers → 50000+ | Longer chains needed for convergence with many parameters |
| `burnIn` | 5000 | Slow convergence in trace plots → 10000+ | Discard pre-convergence samples |
| `thin` | 10 | Very long chains → 20-50; short chains → 5 | Reduces autocorrelation, saves memory |
| `pi` (BayesB/C) | 0.95 | Few known large QTL → 0.99; highly polygenic → 0.90 | \(\pi\) = prior probability a marker has NO effect |
| `df0` (BayesA/B) | 5 | Informative prior on degree-of-freedom → adjust | Higher df0 = prior closer to normal distribution |
| `model` | BayesB | Highly polygenic trait → BayesA; unstable MCMC → BayesC | BayesB can mix poorly with many small-effect QTL |

---

## When Bayes > GBLUP

| Scenario | Recommended | Rationale |
|----------|:---:|------|
| Marker density >50K | **BayesB/C** | Variable selection leverages LD with causal variants |
| Trait with known major QTL | **BayesB** | Explicitly models zero-effect markers |
| n > 500, markers >100K | **BayesC** | More stable MCMC than BayesB in high dimensions |
| Small n (<200) | GBLUP/rrBLUP | Not enough data for variable selection to work reliably |
| Low-density markers | GBLUP (rrBLUP) | Bayes has no advantage when markers don't tag QTL well |
| Published benchmark (plant breeding) | — | BayesB 2-5% accuracy gain over GBLUP at high density |

**Computational note:** BGLR runtime scales roughly linearly with nIter x n_markers. For a 100K-marker panel with 500 individuals, 20000 iterations takes ~30-60 minutes on a modern workstation.

---

## Plant-Specific Notes

### Inbred (self-pollinated) crops — rice, wheat, soybean, barley
- **Prior choice:** BayesC often preferred over BayesB for inbreds
  - Most QTL are already fixed within elite germplasm, so the "spike-and-slab" of BayesB can be too sparse
  - BayesC's common variance assumption is more stable
- **Accuracy expectations:** Gains over GBLUP are modest (2-3%) for elite-by-elite crosses where QTL are largely fixed
- **Model comparison:** Always run rrBLUP baseline first; only invest in Bayes if the trait shows oligogenic architecture

### Outcrossing species — maize, sunflower, brassica, fruit trees
- **Prior choice:** BayesB typically preferred
  - Higher LD decay means more markers needed to tag causal variants
  - Variable selection can identify markers in LD with major QTL
- **Dominance modeling:** Can add dominance ETA term for hybrid prediction:
  ```r
  ETA <- list(
    list(X = geno_add, model = "BayesB"),     # additive
    list(X = geno_dom, model = "BayesB")       # dominance
  )
  ```
- **Heterotic groups:** Run separate within-group predictions; accuracy drops in across-group prediction

### Polyploid species — potato, wheat, strawberry
- BGLR handles continuous marker dosages natively — provide 0-4 dosage scores
- Parameter tuning (pi) may need adjustment: polyploids have more alleles per locus, diluting individual marker effects
- Slightly lower pi (0.90-0.95) recommended to allow more markers into the model

### MCMC convergence checks
- **Always** inspect trace plots for varE (residual variance) — it should stabilize
- Run multiple chains (2-3) with different starting values
- For large genomes (>200K markers): consider chromosome-by-chromosome approach

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| MCMC fails to converge (varE still trending at end) | burnIn too short, nIter too small | Increase both; check trace plots before trusting results |
| All marker effects near zero | pi too high (too sparse); or genuine no-association signal | Lower pi from 0.99 to 0.95; verify trait heritability |
| BGLR hangs or runs extremely slowly | >200K markers without thinning consideration | Use marker pre-selection (e.g., top 50K by GWAS p-value); or increase thinning |
| "Error in solve.default()" | Singular genotype matrix (e.g., identical columns) | Remove duplicate markers; LD-prune to unique sites |
| GEBV accuracy worse than GBLUP | Too few markers for variable selection; overfitting | Fall back to rrBLUP/GBLUP; report honestly |
| Memory allocation error | Genotype matrix too large for RAM | Use `bigmemory` package; or chromosome-by-chromosome BGLR |
| Negative genetic variance estimate | MCMC not converged; model mis-specified | Longer burn-in; switch from BayesB to BayesC |
