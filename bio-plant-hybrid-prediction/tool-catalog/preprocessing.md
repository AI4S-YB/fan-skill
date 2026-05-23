# Hybrid Prediction Preprocessing

**Goal:** Quality control for genotype and phenotype data used in hybrid prediction
**Best for:** Preparing plant breeding data for GCA/SCA analysis and genomic prediction

## Prerequisites

- R 4.0+ with snpReady, data.table, ggplot2
- Genotype data (SNP matrix or VCF)
- Phenotype data (field trial results)
- Pedigree information

## Genotype QC

```r
library(snpReady)

# Load genotype data
M <- as.matrix(read.csv("genotype_matrix.csv", row.names = 1))

# Step 1: Filter markers by call rate
marker_call_rate <- colMeans(!is.na(M))
markers_keep <- marker_call_rate >= 0.90
M <- M[, markers_keep]
print(paste("Markers after call rate filter:", sum(markers_keep)))

# Step 2: Filter individuals by call rate
ind_call_rate <- rowMeans(!is.na(M))
inds_keep <- ind_call_rate >= 0.90
M <- M[inds_keep, ]
print(paste("Individuals after call rate filter:", sum(inds_keep)))

# Step 3: Filter by MAF
p <- colMeans(M, na.rm = TRUE) / 2  # Allele frequency
maf <- pmin(p, 1 - p)
markers_maf <- maf >= 0.05
M <- M[, markers_maf]
print(paste("Markers after MAF filter:", sum(markers_maf)))

# Step 4: LD pruning (optional, for heterotic grouping)
library(bigsnpr)
ld <- snp_ldsc(M)
# Keep markers with r2 < 0.8 in sliding window

# Step 5: Impute missing genotypes
M_imputed <- apply(M, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})
```

## Phenotype QC

```r
# Step 1: Check phenotype distribution
library(ggplot2)
ggplot(hybrid_data, aes(x = Yield)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = mean(hybrid_data$Yield, na.rm = TRUE),
             color = "red", linetype = "dashed") +
  labs(title = "Phenotype Distribution", x = "Yield") +
  theme_minimal()

# Step 2: Remove outliers (>3 SD from mean per environment)
hybrid_data_clean <- hybrid_data %>%
  group_by(Env) %>%
  mutate(
    z_score = (Yield - mean(Yield, na.rm = TRUE)) / sd(Yield, na.rm = TRUE)
  ) %>%
  filter(abs(z_score) <= 3) %>%
  ungroup()

# Step 3: Check for balanced design
with(hybrid_data_clean, table(Male, Female))

# Step 4: Heritability check
library(sommer)
h2_model <- mmer(
  fixed = Yield ~ Env,
  random = ~ Male + Female + Male:Female,
  rcov = ~ units,
  data = hybrid_data_clean
)

# Extract variance components
var_comp <- summary(h2_model)$varcomp
# H2 = Vg / (Vg + Ve / n_rep)
```

## Design Matrix Construction

```r
# Build design matrices for NCII or diallel

# NCII: Males from group A, Females from group B
build_nc2_design <- function(data) {
  data$Male <- factor(data$Male)
  data$Female <- factor(data$Female)

  Z_male <- model.matrix(~ 0 + Male, data = data)
  Z_female <- model.matrix(~ 0 + Female, data = data)
  Z_cross <- model.matrix(~ 0 + Male:Female, data = data)

  return(list(
    Z_male = Z_male,
    Z_female = Z_female,
    Z_cross = Z_cross
  ))
}
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| Marker call rate | > 0.90 | Remove poorly genotyped markers |
| Individual call rate | > 0.90 | Remove poorly genotyped lines |
| MAF threshold | > 0.05 | Remove rare variants |
| Phenotype SD filter | |z-score| <= 3 | Remove extreme outliers |
| Min observations per parent | >= 3 crosses | Stable GCA estimates |
| Min cross per environment | >= 2 replicates | Reliable phenotype |

## Plant-Specific Notes

- Multi-environment trials: check for spatial heterogeneity (row/column effects)
- Multi-year data: year effect is large in most crops — model as fixed or random
- Check for trial management effects (irrigation blocks, fertilizer zones)
- For seed traits: check if seed source (production environment) affects phenotype
- Remove crosses with obvious disease/abiotic damage noted in field books

## Quality Report Generation

```r
generate_qc_report <- function(M, pheno) {
  report <- list()

  # Genotype summary
  report$n_individuals <- nrow(M)
  report$n_markers <- ncol(M)
  report$mean_call_rate <- mean(!is.na(M))
  report$median_maf <- median(pmin(colMeans(M, na.rm = TRUE)/2,
                                   1 - colMeans(M, na.rm = TRUE)/2))

  # Phenotype summary
  report$n_crosses <- nrow(pheno)
  report$n_environments <- length(unique(pheno$Env))
  report$n_males <- length(unique(pheno$Male))
  report$n_females <- length(unique(pheno$Female))
  report$crosses_per_male <- mean(table(pheno$Male))
  report$crosses_per_female <- mean(table(pheno$Female))

  # Print report
  for (key in names(report)) {
    cat(sprintf("%-25s: %s\n", key, report[[key]]))
  }
}
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Genotype and phenotype IDs mismatch | Different naming in datasets | Standardize IDs, check aliases |
| Too few markers after QC | Low-quality genotyping | Relax call rate or use imputation |
| Phenotype data highly unbalanced | Incomplete trial data | Use REML which handles imbalance |
| Parent with zero crosses in design | Incomplete phenotype file | Verify phenotype IDs match genotype IDs |
