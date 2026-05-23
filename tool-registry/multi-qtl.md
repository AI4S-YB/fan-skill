# Multiple QTL Mapping — stepwiseqtl, fitqtl, refineqtl & Multi-Environment QTL

**Goal:** Build multi-QTL models automatically, estimate QTL effects, and perform multi-environment joint QTL analysis
**Best for:** Populations with multiple QTL expected, high-density markers, and multi-environment data

## Prerequisites
- Genetic map built and genotype probabilities calculated
- Phenotype data (single or multi-environment)
- Initial QTL scan completed (for stepwiseqtl starting model)

## Method 1: Stepwise QTL Selection (MQM in R/qtl)

```r
library(qtl)

# Calculate genotype probabilities
mycross <- calc.genoprob(mycross, step = 1, error.prob = 0.001)

# Run stepwise QTL selection
# This automatically adds/drops QTL using forward/backward selection
out.sqtl <- stepwiseqtl(
  cross = mycross,
  pheno.col = "trait1",
  method = "hk",
  model = "normal",
  max.qtl = 10,               # Maximum QTL to consider
  additive.only = TRUE,       # TRUE for RIL/DH, FALSE for F2 (dominance)
  scan.pairs = FALSE,         # Set TRUE to scan for epistatic pairs
  penalties = calc.penalties(operm, alpha = 0.05)  # Penalty thresholds
)

# View selected QTL
summary(out.sqtl)

# Get the QTL object
qtls <- attr(out.sqtl, "qtl")
```

### Penalty Calculation

```r
# Calculate penalties from permutation results
operm <- scanone(mycross, method = "hk", pheno.col = "trait1",
                 n.perm = 1000, verbose = TRUE)
penalties <- calc.penalties(operm, alpha = 0.05)

# penalties contains:
#   main: penalty for adding a main-effect QTL
#   heavy: penalty for adding an epistatic pair
#   light: penalty for re-evaluating a previously dropped QTL
```

## Method 2: Manual Multiple QTL Fitting

```r
# After identifying QTL positions from CIM or stepwiseqtl
qtls <- makeqtl(mycross, chr = c("1", "3", "5", "7"),
                pos = c(45.5, 78.2, 12.8, 92.1),
                what = "prob")

# Fit full model
fit <- fitqtl(
  cross = mycross,
  pheno.col = "trait1",
  qtl = qtls,
  method = "hk",
  model = "normal",
  get.ests = TRUE,    # Get effect estimates
  dropone = TRUE       # Test significance of each QTL in full model
)

# Full model results
summary(fit)

# Interpretation:
# - "%var": total phenotypic variance explained by all QTL
# - "dropone LOD": significance of each QTL conditional on others
# - "ests": additive effect (and dominance for F2)
# - "SEs": standard errors of effect estimates
```

## Method 3: QTL Position Refinement

```r
# Refine QTL positions (iterate to find optimal positions)
refined <- refineqtl(
  cross = mycross,
  pheno.col = "trait1",
  qtl = qtls,
  method = "hk",
  model = "normal",
  keeplodprofile = TRUE,  # Keep LOD profiles for plotting
  verbose = TRUE
)

# Plot LOD profiles for each QTL
plotLodProfile(refined, main = "QTL LOD Profiles (Multiple QTL Model)")

# Updated QTL positions
refined$pos

# Fit final refined model
fit_refined <- fitqtl(mycross, pheno.col = "trait1",
                      qtl = refined, method = "hk",
                      get.ests = TRUE, dropone = TRUE)
summary(fit_refined)
```

## Method 4: Bayesian Interval Mapping (BIM)

```r
# Bayesian QTL mapping — samples from posterior QTL model space
out.bim <- stepwiseqtl(
  cross = mycross,
  pheno.col = "trait1",
  method = "imp",         # Imputation method
  model = "normal",
  max.qtl = 6,
  keeplodprofile = TRUE,
  keeptrace = TRUE        # Keep MCMC trace for diagnostics
)

# Posterior probability of QTL on each chromosome
plot(out.bim)

# MCMC diagnostic plots
plot(attr(out.bim, "trace"))
```

## Multi-Environment QTL Analysis

### R/qtl Approach: Joint Analysis

```r
# For multi-environment data: phenotype matrix with columns per environment
# pheno_cols = c("trait_E1", "trait_E2", "trait_E3")

# Option A: Run per environment and compare
qtl_list <- list()
for (env in c("trait_E1", "trait_E2", "trait_E3")) {
  out <- cim(mycross, pheno.col = env, method = "em",
             n.marcovar = 5, window = 10)
  qtl_list[[env]] <- summary(out, threshold = 2.5)
}

# Find overlapping QTL across environments
# QTL detected in >= 2 environments = "stable QTL"

# Option B: Joint model (requires custom scripting)
# Use mean across environments or fit multivariate mixed model
y_mean <- rowMeans(pull.pheno(mycross)[, c("trait_E1", "trait_E2", "trait_E3")])
mycross$pheno$trait_mean <- y_mean

# Run QTL on mean phenotype
out.mean <- cim(mycross, pheno.col = "trait_mean",
                method = "em", n.marcovar = 5, window = 10)
```

### IciMapping Approach: MET Module

For multi-environment data with QTLxE estimation, IciMapping's MET module is the recommended tool. See `tool-catalog/qtl-icimapping.md` for details.

Advantages of IciMapping MET over per-environment R/qtl:
- Estimates QTLxE interaction variance directly
- Joint model has higher statistical power
- Identifies environment-specific QTL automatically

## Plant-Specific Notes

- **MQM vs CIM**: MQM (stepwiseqtl) typically finds 10-20% more QTL than CIM in plant studies, but at the cost of a slight increase in false positives. Use MQM for exploratory analysis and CIM for confirmatory analysis.
- **Penalty thresholds**: The default penalties from `calc.penalties()` can be conservative (too strict) for plant populations with high LD. Consider using alpha = 0.10 for suggestive QTL discovery.
- **Multi-trait QTL**: If you have multiple correlated traits (e.g., yield and biomass), consider joint multi-trait analysis. In R/qtl, use `scanone` with `phe = 1:n_traits` for multi-trait permutation thresholds.
- **QTLxE decomposition**: In plant breeding, the proportion of QTLxE variance relative to QTL main effect informs selection strategy:
  - QTLxE < QTL main → QTL is stable, good for broad adaptation
  - QTLxE > QTL main → QTL is environment-specific, use for targeted breeding

## Reporting QTL Results

Standard format for plant QTL papers:

| QTL Name | Chr | Peak (cM) | CI (cM) | LOD | PVE (%) | Add | Env-Specific |
|----------|-----|-----------|---------|-----|---------|-----|-------------|
| qPH-1-1 | 1 | 45.5 | 40-51 | 5.2 | 12.3 | -2.8 | -- |
| qPH-3-1 | 3 | 78.2 | 72-85 | 4.1 | 8.9 | 1.5 | -- |
| qPH-5-1 | 5 | 12.8 | 8-18 | 4.8 | 10.1 | -1.9 | E2 only |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `"stepwiseqtl: no QTL found"` | Permutation penalties too strict | Lower alpha or specify manual penalties |
| `"fitqtl: model matrix singular"` | Collinear QTL (same position) | Remove duplicate QTL; check for identical marker positions |
| `"refineqtl: position did not change"` | Already at the peak | Normal behavior; QTL may already be at optimal position |
| Drop-one LOD negative | QTL explains no variance when others in model | QTL may be spurious or fully confounded |
| QTLxE interaction negative variance | Model misspecification | Check phenotype scaling; consider fixed-effect model |
