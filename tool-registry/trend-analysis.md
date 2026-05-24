# Trend Analysis: STEM and Mfuzz-Based Pattern Detection

**Goal:** Identify statistically significant temporal expression trends
(monotonic up, monotonic down, transient, biphasic) and assign genes
to trend profiles.

**Best for:** Dense time courses (>= 6 time points for STEM) where
identifying specific temporal shapes (not just similarity) is the goal.

## Prerequisites

- R 4.0+, `STEM` for STEM analysis (via Java) or `Mfuzz` for profile clustering
- Smoothed, standardized expression matrix
- Time points as ordered numeric vector

## STEM: Short Time-series Expression Miner

STEM is ideal for >= 6 time points and tests whether each predefined
temporal profile is statistically overrepresented.

```r
# STEM requires data in specific format
# Write expression data for STEM (Java-based, runs via command line)

# Prepare input file
stem_input <- data.frame(
  Gene = rownames(z_mat),
  z_mat,
  stringsAsFactors = FALSE
)
write.table(stem_input, "stem_input.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Run STEM from command line:
# java -mx1024M -jar STEM.jar -b stem_input.txt -m 1 -t 0

# --- R-based alternative using Mfuzz profiles ---
library(Mfuzz)

# 1. Define model profiles (normalized time pattern templates)
n_tp <- length(time_points)
model_profiles <- rbind(
  rep(c(1, -1), length.out = n_tp),          # oscillating
  seq(-1, 1, length.out = n_tp),             # monotonic up
  seq(1, -1, length.out = n_tp),             # monotonic down
  c(seq(-1, 1, length.out = n_tp/2),
    seq(1, -1, length.out = n_tp/2)),        # up then down (transient)
  c(seq(1, -1, length.out = n_tp/2),
    seq(-1, 1, length.out = n_tp/2))         # down then up (biphasic)
)

# 2. Correlate each gene with each model profile
cor_matrix <- cor(t(z_mat), t(model_profiles))

# 3. Assign gene to best-matching profile
best_profile <- apply(cor_matrix, 1, which.max)
best_cor     <- apply(cor_matrix, 1, max)

# 4. Count genes per profile
profile_counts <- table(best_profile)
print(profile_counts)

# 5. Permutation test for profile enrichment
n_perm <- 1000
perm_counts <- matrix(0, nrow = n_perm, ncol = nrow(model_profiles))
set.seed(42)
for (p in 1:n_perm) {
  perm_mat <- t(apply(z_mat, 1, sample))  # shuffle time points per gene
  perm_cor <- cor(t(perm_mat), t(model_profiles))
  perm_best <- apply(perm_cor, 1, which.max)
  perm_counts[p, ] <- tabulate(perm_best, nbins = nrow(model_profiles))
}

# Empirical p-value per profile
empirical_p <- sapply(1:nrow(model_profiles), function(i) {
  mean(perm_counts[, i] >= profile_counts[i])
})
names(empirical_p) <- c("oscillating", "monotonic_up", "monotonic_down",
                         "transient_up_down", "biphasic_down_up")
print(round(empirical_p, 4))

# 6. Significant profiles (FDR-corrected)
sig_profiles <- which(p.adjust(empirical_p, "BH") < 0.05)
cat("Significantly enriched profiles:", names(empirical_p)[sig_profiles], "\n")
```

## Mfuzz Trend Detection (< 6 time points)

For fewer time points, use Mfuzz cluster centers as trend profiles
and assign genes to the closest center.

```r
library(Mfuzz)

# After Mfuzz clustering (see mfuzz.md for full pipeline):
# cl <- mfuzz(eset, centers = m1, m = 1.5)
centers <- cl$centers

# Label each cluster center by its trend pattern
# Simple heuristic based on linear regression slope
trend_labels <- apply(centers, 1, function(row) {
  fit <- lm(row ~ time_points)
  slope <- coef(fit)[2]
  r2 <- summary(fit)$r.squared

  if (r2 < 0.3) return("complex")
  if (slope > 0.02) return("up")
  if (slope < -0.02) return("down")
  return("flat")
})

data.frame(cluster = 1:nrow(centers),
           trend = trend_labels,
           n_genes = table(primary_cluster))
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| n_profiles | 50 (STEM) | Number of model profiles |
| correlation_method | pearson | Correlation method for profile matching |
| n_permutations | 1000 | Permutations for significance testing |
| FDR_threshold | 0.05 | Significance cutoff for enriched profiles |
| slope_threshold | 0.02 | Slope cutoff for up/down trend labeling |

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| n_profiles | 50 | If no profiles enriched, reduce to 20-30; if all genes concentrated in one profile, increase to 80-100 | Too many profiles dilutes enrichment power; too few merges biologically distinct patterns |
| n_permutations | 1000 | Increase to 5000-10000 for publication figures | Higher permutation counts stabilize empirical p-values, especially for small gene sets |
| FDR_threshold | 0.05 | Relax to 0.1 for exploratory analysis of noisy plant stress data | Plant stress responses have high biological variability; strict FDR may discard real but subtle trends |
| slope_threshold | 0.02 | Raise to 0.05 for noisy datasets (field-grown plants); lower to 0.01 for tightly controlled growth chamber experiments | Prevents minor expression fluctuations from being labeled as biologically meaningful trends |
| correlation_method | pearson | Switch to spearman for non-linear expression trajectories | Spearman captures monotonic trends that pearson misses (common in developmental gradients) |

## Plant-Specific Notes

- **Circadian oscillating profiles**: add a sine-wave template
  `sin(2*pi*time_points/24 + phase)` to model_profiles for circadian detection.
- **Developmental gradients**: monotonic up/down profiles are expected
  and should be the baseline, not a "finding." Focus on transient
  and biphasic profiles as they indicate developmental transitions.
- **Stress time courses**: rapid-up-then-down (transient) profiles
  are characteristic of early-response genes (e.g., transcription factors).
  Sustained-up profiles indicate adaptation genes.
- **Leaf senescence**: monotonic-down in photosynthetic genes plus
  monotonic-up in proteolysis/autophagy genes is the canonical
  senescence signature.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| No significantly enriched profiles | Sample size too small or data too noisy | Reduce n_profiles (more genes per profile); increase smoothing |
| All genes match one profile | Profile set too coarse or skewed | Add more diverse model profiles |
| STEM Java heap space error | Input matrix too large | Increase -mx memory; or filter to top variable genes |
| Profile correlation matrix all NaN | Zero-variance time points | Check for time points with identical expression across all genes |
| Circular profile assignment | Genes match multiple profiles equally | Use highest absolute correlation; flag ambiguous genes |
