# G×E Modeling (基因型×环境互作建模)

**Goal:** Model genotype-by-environment interaction to understand G×E patterns
and predict genotype performance across environments

## Prerequisites
- R 4.0+, packages: `asreml` or `sommer` or `lme4`, `metan`, `agricolae`
- Multi-environment phenotype data (trait values per genotype × environment)
- Balanced or near-balanced data preferred

## Method 1: Factor-Analytic (FA) Model — Best for ≥5 Env, ≥50 Genotypes

### FA model with ASReml-R

```r
library(asreml)

# Phenotype data: columns = Genotype, Environment, Rep, Yield
data <- read.csv("met_phenotype.csv")

# FA(1) model
fa1 <- asreml(
  fixed = Yield ~ Environment,
  random = ~ fa(Environment, 1):Genotype,
  residual = ~ dsum(~ units | Environment),
  data = data
)

# FA(2) model (most common)
fa2 <- asreml(
  fixed = Yield ~ Environment,
  random = ~ fa(Environment, 2):Genotype,
  residual = ~ dsum(~ units | Environment),
  data = data
)

# Compare models with AIC/BIC
summary(fa1)$aic
summary(fa2)$aic

# Extract genetic correlations between environments
fa2_summary <- summary(fa2)
# The factor loadings represent environmental grouping patterns

# Predict BLUPs
blups <- predict(fa2, classify = "Environment:Genotype")$pvals
```

### FA model with Sommer (open-source alternative)

```r
library(sommer)

# FA(1) with sommer
fa1_sommer <- mmer(
  fixed = Yield ~ Environment,
  random = ~ vsr(Genotype, Gu = fa(Environment, 1)),
  rcov = ~ vsr(units, Gti = diag(Environment)),
  data = data
)
```

## Method 2: AMMI/GGE Biplot — Best for ≥2 Environments

```r
library(metan)

# AMMI analysis
ammi_model <- perf_ammi(
  data,
  env = Environment,
  gen = Genotype,
  rep = Rep,
  resp = Yield,
  verbose = FALSE
)

# AMMI table
ammi_model$ANOVA

# GGE model
gge_model <- gge(
  data,
  env = Environment,
  gen = Genotype,
  resp = Yield
)

# Plot GGE biplot
plot(gge_model, type = 1)  # Basic biplot
plot(gge_model, type = 2)  # Mean vs. Stability
plot(gge_model, type = 3)  # Which-Won-Where
plot(gge_model, type = 4)  # Discriminativeness vs. Representativeness
plot(gge_model, type = 8)  # Ranking genotypes

# Extract GGE scores
gge_scores <- gge_model$gen_effects
```

### AMMI-based stability indices

```r
library(metan)

# Multiple AMMI stability indices
ammi_stab <- ammi_indexes(ammi_model)

# Includes:
# - ASV (AMMI Stability Value)
# - WAASB (Weighted Average of Absolute Scores)
# - EV (Eigenvalue)
# - SIPC (Sum of IPC scores)
```

## Method 3: Reaction Norm — When Environmental Index Is Available

```r
library(lme4)
library(ggplot2)

# Load phenotype data and environmental index
pheno <- read.csv("met_phenotype.csv")
env_idx <- read.csv("env_index.csv")

# Merge
dat <- merge(pheno, env_idx, by = "Environment")

# Fit random regression (reaction norm) model
# Random intercept + random slope on PC1
rn_model <- lmer(
  Yield ~ PC1 + (1 + PC1 | Genotype),
  data = dat
)

summary(rn_model)

# Extract genotype-specific slopes (plasticity)
genotype_slopes <- ranef(rn_model)$Genotype
genotype_slopes$Plasticity <- genotype_slopes$PC1
genotype_slopes$Genotype <- rownames(genotype_slopes)

# Plot reaction norms
ggplot(dat, aes(x = PC1, y = Yield, color = Genotype)) +
  geom_point(alpha = 0.3) +
  stat_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  labs(x = "Environmental Index (PC1)", y = "Yield") +
  theme_minimal()
```

## Method Selection Guide

| Scenario | Recommended | Alternative |
|----------|-------------|-------------|
| ≥5 env, ≥50 gen, need BLUP | FA(2) with ASReml | Sommer |
| 2-4 env, visualization focus | GGE biplot | AMMI |
| Environmental index exists, continuous gradient interest | Reaction Norm | FA + factor regression |
| Unbalanced data, many missing | Mixed model (UNST) | FA (may not converge) |
| Budget-friendly (no ASReml license) | Sommer or metan | lme4 for simple cases |

## Plant Relevance

- **Heritability varies by environment**: Always compute per-environment
  heritability. FA models capture this through heterogeneous residual
  variances.
- **Connectedness is critical**: If certain genotypes appear in only one
  environment, the genetic correlation between that environment and others
  cannot be estimated. Check the connectivity matrix before analysis.
- **Cross-validation**: For breeding applications, validate FA model
  predictions via cross-validation (leave-one-environment-out or
  leave-one-genotype-out).
- **G×E significance is NOT enough**: Always decompose G×E further.
  A significant G×E tells you interaction exists; it doesn't tell you
  which genotypes are interacting with which environments or why.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| FA model fails to converge | Too few environments or genotypes | Reduce FA order or switch to AMMI |
| Negative variance estimates | Sampling variance or model misspecification | Use `asreml` with `!GU` constraint or Bayesian approach |
| AMMI axes not significant | G×E is weak or noise-dominated | Only interpret significant axes (use F-test via `metan`) |
| Reaction norm slope near 0 | PC1 does not capture relevant environmental gradient | Check PC1 biological meaning; try PC2 |
