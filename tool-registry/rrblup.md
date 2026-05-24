# rrBLUP — Ridge Regression Genomic Prediction

**Goal:** Genomic prediction (GS) using ridge regression BLUP — the simplest and most robust GS model
**Approach:** Build additive genomic relationship matrix (G) → mixed.solve with G as random effect covariance → extract GEBVs and assess accuracy
**Best for:** Small breeding populations (<200), all marker densities, baseline accuracy before complex models

## Prerequisites
- R 4.0+, rrBLUP package
- Numeric genotype matrix (samples x markers, coded 0/1/2 or -1/0/1)
- Phenotype vector with matching sample order
- For cross-validation: populate package or custom CV fold assignments

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(rrBLUP)

# ── Step 1: Load and prepare data ──
# geno: samples × markers, coded 0/1/2 (0 = homozygous ref, 1 = het, 2 = homozygous alt)
# pheno: named numeric vector, names must match rownames(geno)
geno <- ${GENO_MATRIX}       # numeric matrix
pheno <- ${PHENO_VECTOR}     # named numeric vector

# ── Step 2: Build additive genomic relationship matrix ──
# Centering: subtract 1 for -1/0/1 coding (standard for rrBLUP)
# For 0/1/2 coding, subtract 1; for already -1/0/1, skip subtraction
G <- A.mat(t(geno) - 1)      # t() because A.mat expects markers × samples

# ── Step 3: Fit mixed model ──
# mixed.solve: y = Xβ + Zu + ε, where u ~ N(0, Kσ²ᵤ)
ans <- mixed.solve(
  y = pheno,
  Z = ${Z_MATRIX},           # design matrix for random effects (default: identity)
  K = G,                     # genomic relationship matrix
  X = ${FIXED_MATRIX},       # fixed effects design matrix (default: intercept only)
  method = "REML"            # REML is standard for variance component estimation
)

# ── Step 4: Extract results ──
gebv <- ans$u                # genomic estimated breeding values (random effects)
gblup_pred <- ans$u + ans$beta  # total predicted value = GEBV + fixed intercept
h2_g <- ans$Vu / (ans$Vu + ans$Ve)  # genomic heritability

# ── Step 5: Accuracy assessment (with independent validation) ──
# When validation set available:
accuracy <- cor(gblup_pred[${VALIDATION_IDS}], pheno[${VALIDATION_IDS}])

# ── Step 6: Cross-validation for accuracy estimation ──
# k-fold CV when no independent validation set
n_folds <- ${N_FOLDS}           # typically 5 or 10
n <- length(pheno)
folds <- sample(rep(1:n_folds, length.out = n))
cv_pred <- numeric(n)

for (k in 1:n_folds) {
  train_idx <- which(folds != k)
  test_idx  <- which(folds == k)
  ans_cv <- mixed.solve(y = pheno[train_idx], K = G[train_idx, train_idx])
  # Predict test set: GEBV for test individuals = G[test, train] %*% alpha
  # where alpha = G⁻¹[test,train] difficult; use kin.blup as alternative:
  ans_cv2 <- kin.blup(
    data = data.frame(id = names(pheno), y = pheno),
    geno = "${GENO_COLUMN}",
    pheno = "${PHENO_COLUMN}",
    K = G
  )
  cv_pred[test_idx] <- ans_cv2$pred[test_idx]
}
cv_accuracy <- cor(cv_pred, pheno)
cat(sprintf("CV accuracy (r): %.3f\n", cv_accuracy))
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `G matrix centering` | `A.mat(t(geno) - 1)` | Already -1/0/1 coded → skip `-1` | rrBLUP expects marker coding centered on zero |
| `method` | REML | Unbalanced design with fixed effects → consider ML | REML gives unbiased variance estimates |
| `N_FOLDS` | 5 | Very small n (<50) → use 10 or leave-one-out | More folds = more training data per fold |
| `Z matrix` | Identity (default) | Multiple records per genotype → custom incidence matrix | Links repeated observations to same genetic effect |
| `X matrix` | Intercept only (default) | Known environmental gradient → include as fixed covariate | Avoids confounding genetic and environmental effects |

---

## When rrBLUP vs GBLUP vs Bayes

| Scenario | Recommended method | Rationale |
|----------|:---:|------|
| n < 200, any marker density | **rrBLUP** | Most robust with limited samples |
| n > 500, moderate markers (10K-50K) | GBLUP | Faster computation, same math |
| High-density markers (50K+), large n | BayesB/BayesC | Variable selection gains 2-5% accuracy |
| Quick feasibility check | **rrBLUP** | Simplest, fastest baseline |
| rrBLUP and GBLUP are equivalent | — | Same accuracy when using identical G matrix |

**Note:** `kin.blup()` is a convenience wrapper that combines A.mat + mixed.solve. Use it for quick analysis; use the two-step approach above when you need to inspect or modify G.

---

## Plant-Specific Notes

### Inbred crops (rice, wheat, soybean, barley, most self-pollinated species)
- Use additive relationship matrix (A.mat) only — **dominance modeling is rarely needed**
- Inbred lines are homozygous, so dominance variance is near zero
- rrBLUP captures additive effects which are the primary target for pure-line selection

### Outcrossing species (maize, sunflower, brassica vegetables, perennial fruit trees)
- Consider dominance matrix (`D.mat()`) if hybrid performance prediction is the goal
- Additive + dominance model: `mixed.solve(y = pheno, K = list(G_add, G_dom))`
- But for line selection: additive-only model is still appropriate

### Polyploid species (potato, wheat, strawberry, sugarcane)
- Coding: 0/1/2 may not capture dosage (tetraploid has 5 possible genotypes)
- Consider allele dosage scores (0/1/2/3/4 for tetraploid) passed as continuous values
- G matrix from dosage scores is approximate but functional for initial assessment

### Low marker density (<1K markers)
- rrBLUP is MORE robust than variable selection methods at low density
- BayesB/C models struggle when causal variants are poorly tagged
- Report accuracy with the caveat that denser markers would improve precision

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "system is computationally singular" | Highly correlated markers, numerical instability in G inverse | LD-prune markers first; add small constant to diagonal: `G <- G + diag(0.01, nrow(G))` |
| "non-conformable arguments" | Genotype and phenotype sample count mismatch | Check `nrow(geno) == length(pheno)`; verify sample IDs match |
| Negative genomic heritability (Vu = 0 or negative) | Low genetic signal relative to error; model mis-specification | Check phenotype distribution; increase population size; try kinship-based model |
| CV accuracy near 0 or negative | Poor marker-trait association; population structure confounded | Check population structure via PCA; use stratified CV folds |
| `A.mat()` returns NaN | Missing genotype values (NA) in matrix | Impute with mean: `geno[is.na(geno)] <- colMeans(geno, na.rm = TRUE)[col(geno)][is.na(geno)]` |
| GEBV all near zero | Shrinkage too strong (small n, weak signal) | Real result — trait may have very low heritability; report honestly |
