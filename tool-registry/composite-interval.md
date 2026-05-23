# Composite Interval Mapping (CIM) — R/qtl cim()

**Goal:** QTL mapping controlling background QTL effects via cofactors
**Best for:** Standard biparental populations (F2, RIL, DH, BC) with moderate-to-high marker density

## Prerequisites
- Genetic map built (est.map or equivalent)
- Genotype probabilities calculated (`calc.genoprob`)
- Phenotype data imported and cleaned

## Basic Usage

```r
library(qtl)

# 1. Calculate genotype probabilities
mycross <- calc.genoprob(mycross, step = 1, error.prob = 0.001)

# 2. Single-QTL scan (IM) — used to select initial cofactors
out.im <- scanone(mycross, method = "em", pheno.col = "trait1")

# 3. Identify cofactors via forward selection
# Option A: Manual cofactor selection
sig_markers <- summary(out.im, threshold = 2.5)  # LOD > 2.5 as cofactor candidates
cofactor_names <- rownames(sig_markers)

# Option B: Automated forward-backward selection
out.cim <- cim(mycross, pheno.col = "trait1",
               method = "em",
               n.marcovar = 5,   # Max number of cofactors
               window = 10)      # Window size (cM) to exclude cofactors near test position

# 4. Permutation test for CIM LOD threshold
operm.cim <- cim(mycross, pheno.col = "trait1",
                 method = "em",
                 n.marcovar = 5,
                 window = 10,
                 n.perm = 1000)

# 5. Identify significant QTL
summary(out.cim, perms = operm.cim, alpha = 0.05, pvalues = TRUE)

# 6. Refine QTL positions
refined <- refineqtl(mycross, pheno.col = "trait1",
                     qtl = makeqtl(mycross, chr = "1", pos = 45.5),
                     method = "hk",
                     model = "normal")
```

## Cofactor Selection Strategy

```r
# Recommended: stepwise forward selection + backward elimination
add_cofactors <- function(cross, pheno.col, n.marcovar = 5) {
  # Forward selection to pick cofactors
  cofactors <- NULL
  for (i in 1:n.marcovar) {
    out <- scanone(cross, pheno.col = pheno.col, method = "em",
                   addcovar = cofactors)
    best <- which.max(out$lod)
    cofactors <- c(cofactors, rownames(out)[best])
  }
  # Backward elimination
  while (length(cofactors) > 0) {
    models <- sapply(seq_along(cofactors), function(j) {
      fitqtl(cross, pheno.col = pheno.col,
             qtl = makeqtl(cross, chr = names(cofactors[-j]),
                           pos = cofactors[-j],
                           what = "prob"),
             method = "hk", get.ests = TRUE,
             dropone = FALSE)$result.full[1, "LOD"]
    })
    # Drop the least important cofactor if LOD drop < threshold
  }
  return(cofactors)
}
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| n.marcovar | 3-5 for <500 markers, 5-10 for >500 | More markers → more potential cofactors |
| window | 10 cM (standard) or 15 cM (sparse map) | Exclude cofactors near the test position |
| method | "em" (accurate) or "hk" (fast, large datasets) | HK recommended for >5000 markers |
| error.prob | 0.001 (SSR), 0.01 (GBS) | Genotyping error rate for genotype probability |

## Model Refinement

After identifying QTL via CIM:

```r
# Fit full QTL model
qtls <- makeqtl(mycross, chr = c("1", "3", "5"),
                pos = c(45.5, 78.2, 12.8),
                what = "prob")

# Fit and get statistics
fit <- fitqtl(mycross, pheno.col = "trait1",
              qtl = qtls, method = "hk",
              model = "normal", get.ests = TRUE,
              dropone = TRUE)

summary(fit)
# Look at:
# - % variance explained per QTL
# - LOD scores
# - Additive and dominance effects
# - Drop-one analysis to see if each QTL is significant in the full model

# Refine positions
refined <- refineqtl(mycross, pheno.col = "trait1",
                     qtl = qtls, method = "hk",
                     model = "normal")

# Plot refined QTL LOD profile
plotLodProfile(refined, main = "Refined QTL LOD Profiles")
```

## CIM vs IM — When to Use Each

| Aspect | IM | CIM |
|--------|----|-----|
| Background QTL control | None | Yes, via cofactors |
| Multiple linked QTL | Can produce ghost QTL | Better discrimination |
| Statistical power | Lower | Higher |
| Model complexity | Simple | Moderate |
| Computation time | Fast | Slower (cofactor selection) |
| Minimum markers | 50+ | 100+ |
| False positive rate | Higher with linked QTL | Better controlled |

## Plant-Specific Notes

- **CIM is the standard for plant QTL studies** — it is the most commonly reported method in crop QTL papers (rice, wheat, maize, soybean).
- **Cofactor window size in plants**: Use 10 cM for most species. Adjust to 15 cM if the map is sparse, or 5 cM for ultra-dense maps (>5000 markers). 
- **Dominance effects in F2**: Use `sim.geno()` instead of `calc.genoprob()` if you need to test for dominance effects via `fitqtl`.
- **Epistasis**: CIM does not model QTLxQTL epistasis. Use `scantwo()` for two-dimensional scans.
- **QTL naming convention**: Plant journals expect QTLs named as qTrait-Chromosome-Number (e.g., qPH-1-1, qPH-1-2 for plant height QTL on chromosome 1).

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `"Need to run calc.genoprob first"` | Missing genotype probabilities | `calc.genoprob()` before cim() |
| `"n.marcovar too large"` | More cofactors than useful markers | Reduce n.marcovar or increase window |
| `"overfitting"` | Too many cofactors for sample size | Rule of thumb: n.marcovar < n_samples / 20 |
| `"different QTL when changing window"` | Window size affecting cofactor exclusion | Test window = 5, 10, 15 and compare results |
| CIM LOD << IM LOD | Cofactors absorbing QTL effect | The QTL may be colocating with cofactor; adjust window |
