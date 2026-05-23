# Plant Hybrid Prediction Specific Considerations

## Breeding Design Types

### NCII (North Carolina Design II)
- **Structure**: Set A males x Set B females (factorial)
- **Typical use**: Rice, maize hybrid breeding
- **Analysis**: GCA + SCA via mixed model (ASReml-R / sommer)
- **Advantage**: Simple design, easy to interpret
- **Disadvantage**: Assumes clear heterotic groups exist

### Diallel (Complete / Half)
- **Structure**: All possible crosses among n parents
- **Griffing's methods**: Method 1 (full diallel with reciprocals), Method 2 (full without reciprocals), Method 3 (half diallel), Method 4 (half diallel with parents)
- **Typical use**: Small-scale combining ability studies (<15 parents)
- **Advantage**: Complete genetic information
- **Disadvantage**: n parents = n(n-1)/2 crosses (scales poorly)

### Single-Cross Evaluation
- **Structure**: Individual F1 crosses from heterotic groups
- **Typical use**: Late-stage hybrid testing
- **Analysis**: BLUP of cross performance across environments
- **Advantage**: Directly evaluates commercial product
- **Disadvantage**: No inference about untested crosses without genomic data

## GCA/SCA Variance Components by Crop

| Crop | Trait | sigma_GCA^2 ratio | sigma_SCA^2 ratio | Interpretation |
|------|-------|-------------------|-------------------|----------------|
| Rice | Grain yield | 0.70-0.85 | 0.15-0.30 | GCA dominates |
| Rice | Grain quality | 0.80-0.95 | 0.05-0.20 | Mostly GCA |
| Maize | Grain yield | 0.45-0.65 | 0.35-0.55 | Both important |
| Maize | Plant height | 0.60-0.80 | 0.20-0.40 | GCA dominant |
| Rapeseed | Seed yield | 0.40-0.60 | 0.40-0.60 | SCA often = GCA |
| Sunflower | Seed yield | 0.50-0.70 | 0.30-0.50 | Mixed |
| Sorghum | Grain yield | 0.55-0.75 | 0.25-0.45 | Slightly GCA dominant |

## Heterotic Group Discovery

### Genomic Methods

```r
# PCA-based grouping
library(snpReady)
G_matrix <- Gmatrix(SNPmatrix = M, method = "VanRaden")
pca <- prcomp(G_matrix)
groups <- kmeans(pca$x[, 1:5], centers = 3)

# Phylogenetic tree
dist_genetic <- as.dist(1 - G_matrix)  # 1 - relationship
tree <- hclust(dist_genetic, method = "ward.D2")
plot(tree)
```

### Validation of Heterotic Groups

The true test of heterotic groups is hybrid performance:
- Within-group crosses: lower heterosis
- Between-group crosses: higher heterosis

```r
# Test: compare between-group vs within-group hybrid performance
hybrid_data$cross_type <- ifelse(
  hybrid_data$male_group != hybrid_data$female_group,
  "between",
  "within"
)
t.test(yield ~ cross_type, data = hybrid_data)
```

## Genomic Prediction Models for Hybrids

### GBLUP with Dominance

```r
library(BGLR)

# Compute additive and dominance relationship matrices
# VanRaden G matrix (additive)
G_A <- (M %*% t(M)) / sum(2 * p * (1 - p))

# Dominance relationship matrix
D <- ifelse(M == 1, -2 * p * (1 - p),
             ifelse(M == 0, -2 * p^2, -2 * (1 - p)^2))
G_D <- (D %*% t(D)) / sum((2 * p * (1 - p))^2)

# Model
# y = mu + GCA_parent1 + GCA_parent2 + SCA_cross + e
```

### Cross-Prediction Accuracy Benchmarks

| Crop | Trait | Prediction accuracy (r) | Model |
|------|-------|------------------------|-------|
| Maize | Grain yield | 0.50-0.70 | GBLUP with dominance |
| Rice | Grain yield | 0.35-0.55 | GBLUP |
| Wheat (hybrid) | Grain yield | 0.40-0.60 | BayesR |
| Rapeseed | Seed yield | 0.30-0.50 | GBLUP |
| Sunflower | Oil yield | 0.45-0.65 | GBLUP |

## Optimal Mate Selection

### OCS with Genetic Algorithm

```r
# Using optiSel package
library(optiSel)

# Calculate optimal contributions
cand <- candes(phen = data, cont = Ped)

# Solve OCS problem
con <- list(uniform = "female", ub = rep(0.1, nrow(cand)))
opt <- opticont(method = "min.IBD", cand = cand, con = con)
```

### Breeding Pipeline Constraints

Common constraints in plant breeding mate selection:
- Maximum number of crosses per parent (e.g., <= 5)
- Minimum number of parents used (e.g., >= 20) — maintain diversity
- Maximum relatedness between selected parents (e.g., coancestry < 0.125)
- Avoid known lethal allele combinations
- Enforce heterotic group structure (A-group male, B-group female)

## Software Ecosystem for Plant Hybrid Prediction

| Software | Language | Best For | Plant Support |
|----------|----------|----------|---------------|
| ASReml-R | R | Mixed models (GCA/SCA), spatial analysis | Yes, plant breeding focused |
| sommer | R | Multivariate mixed models, GxE | Yes, actively developed for plants |
| BGLR | R | Bayesian genomic prediction | Yes, flexible |
| rrBLUP | R | Ridge regression BLUP | Fast, large datasets |
| AlphaPlant | GUI | Complete plant breeding pipeline | Commercial (maize, rice, wheat) |
| BreedBase | Web | Breeding program management | Cassava, sweet potato, rice |
