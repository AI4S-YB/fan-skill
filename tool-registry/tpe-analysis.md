# TPE Analysis (目标环境群体分析)

**Goal:** Delineate Target Population of Environments (TPE) for breeding
strategy optimization

**Best for:** ≥5 environments with environmental covariate data

## Prerequisites
- R 4.0+, packages: `cluster`, `factoextra`, `NbClust`, `tidyverse`, `pheatmap`
- Environmental covariate or genetic correlation matrix for environments

## Basic Usage

### Method 1: Climate-based Clustering

```r
library(cluster)
library(factoextra)
library(tidyverse)

# Load environmental covariates
env_data <- read.csv("env_covariates.csv", row.names = 1)
env_scaled <- scale(env_data)

# Determine optimal K using multiple indices
library(NbClust)
nb <- NbClust(env_scaled, distance = "euclidean",
              min.nc = 2, max.nc = min(8, nrow(env_scaled) - 1),
              method = "kmeans", index = "silhouette")
optimal_k <- nb$Best.nc[1, "Number_clusters"]

# K-means clustering
set.seed(42)
km <- kmeans(env_scaled, centers = optimal_k, nstart = 25)

# Assign TPE labels
env_data$TPE <- paste0("TPE_", km$cluster)

# Silhouette validation
sil <- silhouette(km$cluster, dist(env_scaled))
cat(sprintf("Average silhouette width: %.3f\n", mean(sil[, 3])))
# > 0.3 = acceptable, > 0.5 = good, > 0.7 = strong

# Visualize clusters in PCA space
fviz_cluster(km, data = env_scaled,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal(),
             main = "TPE Delineation by Climate Variables")
```

### Method 2: Genetic Correlation-based Clustering

```r
# Use FA model output: genetic correlation matrix between environments
# gen_cor_mat: environments x environments matrix of genetic correlations

library(pheatmap)

# Compute distance from genetic correlation
gen_dist <- as.dist(1 - gen_cor_mat)  # cor=1 means distance=0

# Hierarchical clustering
hc <- hclust(gen_dist, method = "ward.D2")

# Plot dendrogram
plot(hc, main = "TPE Delineation by Genetic Correlation",
     xlab = "", sub = "")

# Cut into groups
tpe_groups <- cutree(hc, k = optimal_k)

# Heatmap of genetic correlations with TPE grouping
pheatmap(gen_cor_mat,
         clustering_method = "ward.D2",
         display_numbers = TRUE,
         main = "Genetic Correlations Between Environments")
```

### Method 3: Combined Approach

```r
# Combine climate distance and genetic correlation distance
climate_dist <- dist(scale(env_data))
genetic_dist <- as.dist(1 - gen_cor_mat)

# Scale both to [0,1]
climate_dist_scaled <- climate_dist / max(climate_dist)
genetic_dist_scaled <- genetic_dist / max(genetic_dist)

# Weighted combination (adjust alpha based on data quality)
alpha <- 0.5  # 0 = pure genetic, 1 = pure climate
combined_dist <- alpha * climate_dist_scaled +
                 (1 - alpha) * genetic_dist_scaled

hc_combined <- hclust(combined_dist, method = "ward.D2")
plot(hc_combined, main = "TPE: Combined Climate + Genetic Distance")
```

## Cluster Validation

```r
# Bootstrap stability of clusters
library(fpc)
boot_res <- clusterboot(env_scaled,
                        B = 100,
                        clustermethod = kmeansCBI,
                        krange = optimal_k,
                        seed = 42)
# Jaccard similarity for each cluster (should be > 0.75)
boot_res$bootmean

# Compare multiple clustering methods
hc_kmeans <- cutree(hclust(dist(env_scaled), method = "ward.D2"), k = optimal_k)
table(km$cluster, hc_kmeans)  # Agreement between methods
```

## Plant Relevance

- **Dynamic TPEs**: TPE boundaries shift with climate change. Always note the
  time period of your environmental data and re-evaluate periodically.
- **Breeding implications**: If two environments belong to different TPEs with
  low genetic correlation, they may need separate breeding programs.
- **TPE size matters**: A TPE with only 1-2 historical environments has
  limited data. Consider merging with adjacent TPEs if they are too small.
- **Within-TPE variation**: Even within a TPE, there is substantial
  environment-to-environment variation. TPEs are useful abstractions, not
  perfect predictions.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| Silhouette < 0.2 | No genuine clustering structure | Don't force clusters; report continuous gradients |
| All env in one TPE | Homogeneous environments | TPE analysis not informative for this dataset |
| TPE assignment unstable | Few environments or redundant variables | Use bootstrap to assess stability |
