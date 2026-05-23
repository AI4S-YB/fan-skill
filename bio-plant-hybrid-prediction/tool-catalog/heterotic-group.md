# Heterotic Group Classification

**Goal:** Classify inbred lines into heterotic groups using genomic, pedigree, or phenotypic data
**Best for:** Establishing heterotic patterns in plant hybrid breeding programs

## Prerequisites

- Genotype data (SNP matrix: individuals x markers)
- Pedigree data (optional)
- R 4.0+ with snpReady, adegenet, ggtree, cluster

## Genomic Heterotic Grouping

### PCA-Based Method

```r
library(snpReady)

# Build genomic relationship matrix
M <- as.matrix(genotype_data[, -1])  # Remove ID column
G <- Gmatrix(SNPmatrix = M, method = "VanRaden",
             missingValue = NA, thresh.missing = 0.1,
             maf = 0.05)

# PCA on G matrix
pca <- prcomp(G, scale. = TRUE, center. = TRUE)
pca_scores <- as.data.frame(pca$x)

# Determine number of groups (elbow method)
wss <- sapply(1:10, function(k) {
  km <- kmeans(pca_scores[, 1:min(10, ncol(pca_scores))], centers = k, nstart = 25)
  km$tot.withinss
})
plot(1:10, wss, type = "b", xlab = "Number of groups (k)", ylab = "WSS")

# Assign groups (choose k based on elbow)
k_optimal <- 3  # Adjust based on plot
set.seed(42)
groups <- kmeans(pca_scores[, 1:5], centers = k_optimal, nstart = 50)

# Add group labels to data
inbred_lines$heterotic_group <- paste0("HG", groups$cluster)
```

### Phylogenetic Tree Method

```r
library(ape)
library(ggtree)

# Calculate genetic distance
dist_mat <- as.dist(1 - G)  # 1 - relationship = genetic distance
diag(dist_mat) <- 0

# Build NJ tree
nj_tree <- nj(dist_mat)

# Plot with group colors
group_colors <- rainbow(k_optimal)[groups$cluster]
names(group_colors) <- rownames(G)

ggtree(nj_tree, layout = "circular") +
  geom_tippoint(aes(color = heterotic_group), size = 3) +
  theme_tree()
```

## Pedigree-Based Grouping

```r
library(nadiv)

# Build A matrix from pedigree
A_mat <- makeA(pedigree)

# Calculate coancestry distance
coancestry_dist <- as.dist(1 - A_mat)

# Hierarchical clustering
hc <- hclust(coancestry_dist, method = "ward.D2")
plot(hc, labels = pedigree$ID, cex = 0.6)

# Cut into groups
groups_ped <- cutree(hc, k = k_optimal)
```

## Validating Heterotic Groups

```r
# If hybrid performance data is available:
# Compare mean hybrid performance within vs between groups

hybrid_data$cross_type <- with(hybrid_data, {
  male_group <- inbred_lines$heterotic_group[match(Male, inbred_lines$ID)]
  female_group <- inbred_lines$heterotic_group[match(Female, inbred_lines$ID)]
  ifelse(male_group != female_group, "Between", "Within")
})

# Statistical test
library(ggplot2)
ggplot(hybrid_data, aes(x = cross_type, y = Yield, fill = cross_type)) +
  geom_boxplot() +
  labs(title = "Hybrid Performance: Between vs Within Heterotic Groups",
       y = "Yield", x = "Cross Type") +
  theme_minimal()

t.test(Yield ~ cross_type, data = hybrid_data)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| MAF filter | > 0.05 | Remove rare alleles (noise in grouping) |
| Missing call rate | < 10% | Remove poorly genotyped lines/markers |
| LD pruning | r2 < 0.8 | Reduce redundancy in markers |
| Number of PCs | 5-10 | Capture major population structure |
| Silhouette score | > 0.3 | Minimum acceptable clustering quality |

## Plant-Specific Notes

- Major crops (maize, rice, sorghum) have well-established heterotic groups — validate known patterns first
- For new crops without known groups: PCA is the most exploratory approach
- Hybrid performance data is the gold standard for validating groups
- In crops with strong population structure (rice indica/japonica), heterotic groups align with subspecies
- For synthetic populations: genomic grouping may not map to clear heterotic patterns

## Silhouette Analysis

```r
library(cluster)

# Compute silhouette score
sil <- silhouette(groups$cluster, dist(pca_scores[, 1:5]))
mean_sil <- mean(sil[, 3])

plot(sil, border = NA,
     main = paste0("Silhouette Plot (Mean = ", round(mean_sil, 3), ")"))

# If mean_sil > 0.5: reasonable grouping
# If mean_sil 0.3-0.5: weak grouping structure
# If mean_sil < 0.3: no clear groups — consider alternative methods
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| All lines in one group | Insufficient genetic diversity | Add more diverse germplasm, check if PCA shows structure |
| Silhouette < 0.2 | No clear genetic structure | Breeding population may be genetically homogeneous |
| Groups don't match known heterotic pattern | Genomic clustering picks up breeding program structure, not heterotic pattern | Include historical heterotic tester data |
| Number of groups unclear | Continuous population structure | Use DAPC (adegenet) instead of k-means |
