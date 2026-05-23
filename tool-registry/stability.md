# Stability Analysis (基因型稳定性分析)

**Goal:** Identify genotypes with consistent performance across environments
using multiple stability metrics for cross-validation

## Prerequisites
- R 4.0+, packages: `metan`, `agricolae`, `stats`, `tidyverse`
- Multi-environment phenotype data: genotype × environment × trait

## Basic Usage

### Multiple Stability Metrics with metan

```r
library(metan)

# Data format: columns = Genotype, Environment, Rep, Yield
data <- read.csv("met_phenotype.csv")

# Comprehensive stability analysis
stability_res <- ge_stats(
  data,
  env = Environment,
  gen = Genotype,
  rep = Rep,
  resp = Yield,
  verbose = TRUE
)

# Returns stability measures:
# - Y: mean yield per genotype
# - CV: coefficient of variation (static stability, lower = more stable)
# - Var: environmental variance (static stability)
# - Shukla: Shukla stability variance (lower = more stable)
# - Wi: Wricke ecovalence (lower = more stable)
# - b (FW regression): Finlay-Wilkinson regression slope
#   b ≈ 1.0 = average response
#   b < 1.0 = above-average stability in poor environments
#   b > 1.0 = responsive to good environments
# - S2di: Eberhart-Russell deviation from regression (lower = more predictable)
# - superiority measure (Pi): distance from best genotype per environment

# View summary
head(stability_res$stats)
```

### Individual Stability Metrics with agricolae

```r
library(agricolae)

# Prepare data
data_wide <- data %>%
  select(Genotype, Environment, Yield) %>%
  pivot_wider(names_from = Environment, values_from = Yield, values_fn = mean)

# Finlay-Wilkinson regression
fw_res <- stability.par(data_wide[, -1],
                        rep = 3,  # number of reps
                        MSerror = mse,
                        cova = TRUE)
# Extract b values and test if b ≠ 1

# AMMI analysis
ammi_res <- AMMI(ENV = data$Environment,
                 GEN = data$Genotype,
                 REP = data$Rep,
                 Y = data$Yield,
                 MSE = mse,
                 console = TRUE)
```

### Cross-Validate Multiple Stability Metrics

```r
# Extract rankings from different metrics
rankings <- stability_res$stats %>%
  select(GEN, Y, CV, Var, Shukla, Wi, FW, S2di, Pi) %>%
  mutate(
    rank_Y = rank(-Y),           # higher yield = better rank
    rank_CV = rank(CV),          # lower CV = better rank
    rank_Var = rank(Var),        # lower variance = better rank
    rank_Shukla = rank(Shukla),  # lower = better rank
    rank_FW = rank(abs(FW - 1)), # closer to 1 = better rank
    rank_S2di = rank(S2di),      # lower = better rank
    rank_Pi = rank(Pi)           # lower = better rank
  ) %>%
  select(GEN, starts_with("rank_"))

# Kendall's W coefficient of concordance
library(irr)
kendall_w <- kendall(rankings[, -1])
cat(sprintf("Kendall W = %.3f (p = %.4f)\n",
            kendall_w$value, kendall_w$p.value))

# If W > 0.6: good agreement → stability assessment is robust
# If W < 0.4: poor agreement → different metrics capture different concepts
```

## Stability Metric Interpretation

| Metric | Low Value | High Value | Ideal |
|--------|-----------|------------|-------|
| CV | Very stable | Unstable | Low |
| Environmental variance (Var) | Very stable | Highly variable | Low |
| Shukla variance | Stable | Unstable | Low |
| Wricke ecovalence (Wi) | Small G×E contribution | Large G×E contribution | Low |
| FW regression slope (b) | Stability in poor env | Responsive to good env | ≈ 1.0 |
| S2di (deviation) | Predictable response | Unpredictable | Low |
| Pi (superiority) | Close to best genotype | Far from best | Low |

## Selecting Stable Genotypes

```r
# Multi-trait stability index (MTSI)
library(metan)

mtsi_res <- mtsi(
  data,
  env = Environment,
  gen = Genotype,
  rep = Rep,
  resp = everything(),
  # Select genotypes in top 15% for stability
  perc = 15
)

# Visualize selection differential
plot(mtsi_res)
```

## Plant Relevance

- **Breeding objective matters**: Choosing "stable" genotypes depends on
  what stability means for your program:
  - **Subsistence farming**: Static stability (low CV, consistent yield even
    in bad years)
  - **Commercial production**: Dynamic stability (b ≈ 1, responds to inputs
    but predictably)
  - **Broad adaptation**: High mean + low Pi + low Shukla
- **Don't confuse low yield with stability**: A genotype can be "stable"
  simply because it yields poorly everywhere. Always co-evaluate mean
  performance and stability. GGE biplots are ideal for this.
- **Crossover interactions**: When genotype rankings change between
  environments (crossover G×E), stability metrics based on variance alone
  may be misleading. Complement with GGE biplot "which-won-where" analysis.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| High stability = low yield | Confounding stability with poor performance | Co-plot mean vs stability (GGE type 2) |
| Kendall W < 0.3 | Different definitions of stability conflict | Report metric-by-metric, don't force consensus |
| FW regression extrapolation | Extreme environments beyond range | Check environmental index range; flag extrapolation |
