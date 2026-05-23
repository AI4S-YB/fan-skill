# GCA/SCA Estimation

**Goal:** Estimate general combining ability (GCA) and specific combining ability (SCA) from hybrid trials
**Best for:** NCII design, diallel analysis, and factorial mating designs in plant breeding

## Prerequisites

- R 4.0+ with sommer, lme4, or ASReml-R
- Phenotype data (yield, quality, or other traits)
- Experimental design data (male, female, block, environment)
- Pedigree data (optional, for BLUP with relationship matrix)

## Basic Usage (sommer)

```r
library(sommer)

# NCII Design Analysis
# Data format: hybrid_data with columns: Male, Female, Block, Env, Yield

# Model: Yield ~ Env + Male + Female + Male:Female
model_nc2 <- mmer(
  fixed = Yield ~ Env,
  random = ~ Male + Female + Male:Female,
  rcov = ~ units,
  data = hybrid_data
)

# Extract variance components
summary(model_nc2)$varcomp

# Extract BLUPs
gca_male_blup <- ranef(model_nc2)$Male
gca_female_blup <- ranef(model_nc2)$Female
sca_blup <- ranef(model_nc2)$`Male:Female`

# Baker's ratio (proportion of GCA in total genetic variance)
var_gca <- gca_male_var + gca_female_var
var_sca <- sca_var
baker_ratio <- var_gca / (var_gca + var_sca)
```

## Diallel Analysis (Griffing Methods)

```r
# Half diallel (Method 4: F1 only, no parents, no reciprocals)
model_diallel <- mmer(
  fixed = Yield ~ 1,
  random = ~ Parent1 + Parent2 + Parent1:Parent2,
  rcov = ~ units,
  data = diallel_data
)

# Full diallel with reciprocals (Method 1)
model_diallel_full <- mmer(
  fixed = Yield ~ 1,
  random = ~ Parent1 + Parent2 + Parent1:Parent2 + Recip,
  rcov = ~ units,
  data = diallel_full_data
)
```

## Multi-Environment Trial (MET)

```r
# NCII across environments
model_nc2_met <- mmer(
  fixed = Yield ~ Env,
  random = ~ Male + Female + Male:Female +
             Male:Env + Female:Env + Male:Female:Env,
  rcov = ~ vsr(dsr(Env), units),
  data = met_data
)

# Extract environment-specific GCA/SCA
summary(model_nc2_met)$varcomp
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| Heritability | h2 > 0.1 | Below 0.1, GCA estimates unreliable |
| Baker's ratio | > 0.5 | GCA should explain majority of genetic variance |
| Min crosses per parent | >= 3 | Lower = GCA estimate unstable |
| Min environments | >= 2 | Single environment = no GxE estimate |

## Plant-Specific Notes

- For rice: GCA typically dominates (Baker's ratio 0.7-0.85) — selecting parents based on GCA is effective
- For maize: SCA often contributes 30-50% of genetic variance — specific combinations matter
- For polyploids (wheat, potato): include ploidy in kinship matrix calculation
- Multi-year data is critical: single-year GCA estimates can be very unstable
- Male vs female effects: in some crops, reciprocal effects matter (maternal effect, cytoplasmic)

## GCA-Based Parent Selection

```r
# Rank parents by GCA
gca_rank <- data.frame(
  Parent = names(gca_male_blup$Yield),
  GCA = gca_male_blup$Yield
)
gca_rank <- gca_rank[order(gca_rank$GCA, decreasing = TRUE), ]

# Select top 20% as elite parents
n_select <- ceiling(nrow(gca_rank) * 0.2)
elite_parents <- gca_rank$Parent[1:n_select]

print(paste("Baker's ratio:", round(baker_ratio, 3)))
print(paste("Selected", n_select, "elite parents"))
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Model not converging | Unbalanced design, too few replications | Check design matrix, remove parents with <2 crosses |
| "Singular" warning in sommer | Incomplete rank in design matrix | Verify at least 2 levels per random effect |
| Negative variance components | Insufficient data for that component | Fix to 0 or remove random term |
| GCA estimates very different between years | Large GxE interaction | Use multi-year BLUP |
