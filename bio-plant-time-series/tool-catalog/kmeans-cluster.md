# k-Means Clustering for Time-Series Expression

**Goal:** Hard-cluster genes into k groups by temporal expression pattern.

**Best for:** Smaller gene sets (< 1000 genes) where interpretation
benefits from discrete cluster assignment.

## Prerequisites

- R 4.0+, base R `kmeans()` or `stats` package
- Smoothed, normalized expression matrix (gene x time-point)
- Gene filtering: remove constant or near-constant genes

## Basic Usage

```r
# 1. Standardize gene-wise (z-score)
gene_means <- rowMeans(smoothed_mat)
gene_sds   <- apply(smoothed_mat, 1, sd)
z_mat <- (smoothed_mat - gene_means) / gene_sds
z_mat <- z_mat[gene_sds > 0, ]  # remove zero-variance genes

# 2. Determine optimal k: elbow method
set.seed(42)
wss <- sapply(2:15, function(k) {
  kmeans(z_mat, centers = k, nstart = 25)$tot.withinss
})

pdf("kmeans_elbow.pdf")
plot(2:15, wss, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster SS",
     main = "Elbow Plot for k Selection")
dev.off()

# 3. Gap statistic (alternative, more rigorous)
library(cluster)
set.seed(42)
gap <- clusGap(z_mat, FUN = kmeans, K.max = 15, B = 100, nstart = 25)
optimal_k <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "Tibs2001SEmax")
cat("Optimal k by gap statistic:", optimal_k, "\n")

# 4. Run k-means
set.seed(42)
k <- optimal_k  # or from elbow
km <- kmeans(z_mat, centers = k, nstart = 50, iter.max = 100)

# 5. Silhouette score for cluster quality
library(cluster)
sil <- silhouette(km$cluster, dist(z_mat))
cat("Average silhouette width:", mean(sil[, 3]), "\n")
# >= 0.5: good; 0.3-0.5: acceptable; < 0.3: poor

# 6. Extract cluster assignments
cluster_assignments <- data.frame(
  gene = rownames(z_mat),
  cluster = km$cluster,
  stringsAsFactors = FALSE
)

# 7. Plot cluster centers
centers <- km$centers
pdf("kmeans_clusters.pdf", width = 12, height = 8)
par(mfrow = c(ceiling(k/3), 3), mar = c(3, 3, 2, 1))
for (i in 1:k) {
  plot(time_points, centers[i, ], type = "b", pch = 19, col = i + 1,
       xlab = "", ylab = "", main = paste("Cluster", i, "(n =", km$size[i], ")"))
  abline(h = 0, lty = 2, col = "grey")
}
dev.off()

# 8. Per-cluster mean expression profile with SE ribbon
pdf("kmeans_profiles.pdf", width = 10, height = 8)
par(mfrow = c(ceiling(k/3), 3), mar = c(3, 3, 2, 1))
for (i in 1:k) {
  cluster_genes <- names(km$cluster[km$cluster == i])
  cluster_mat <- z_mat[cluster_genes, , drop = FALSE]
  cluster_mean <- colMeans(cluster_mat)
  cluster_se <- apply(cluster_mat, 2, sd) / sqrt(nrow(cluster_mat))

  ylim <- range(c(cluster_mean + cluster_se, cluster_mean - cluster_se))
  plot(time_points, cluster_mean, type = "b", pch = 19, col = i + 1,
       ylim = ylim, xlab = "", ylab = "",
       main = paste("Cluster", i, "(n =", length(cluster_genes), ")"))
  arrows(time_points, cluster_mean - cluster_se,
         time_points, cluster_mean + cluster_se,
         angle = 90, code = 3, length = 0.05, col = i + 1)
  abline(h = 0, lty = 2, col = "grey")
}
dev.off()
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| centers (k) | Gap statistic | Number of clusters |
| nstart | 50 | Random initializations; higher = more stable |
| iter.max | 100 | Maximum iterations per run |
| algorithm | Hartigan-Wong | k-means algorithm variant |
| scaling | z-score | Standardization method before clustering |

## Plant-Specific Notes

- **Small focused gene sets** (e.g., TFs only, hormone pathway genes):
  k = 3-6 is usually optimal regardless of gene count. The biological
  signal in these sets is structured into few patterns.
- **DEG subsets**: cluster differentially expressed genes separately
  from the full transcriptome. k-means on DEGs reveals sub-patterns
  (e.g., early vs late responders) that are diluted in whole-genome clustering.
- **Time-point ordering**: k-means ignores temporal order (unlike spline
  or STEM). If temporal adjacency matters strongly, verify that adjacent
  time points fall in the same half of the cluster center profile.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| `kmeans: empty cluster` | k too large for data | Reduce k or increase nstart |
| `gap statistic fails` | Too few genes (< 50) | Use elbow method or fixed k |
| Silhouette < 0.2 | No real cluster structure | Filter to most variable genes; reduce k |
| All genes in one cluster | k too small or data not clustered | Check variance distribution; try PCA first |
| Centers unstable between runs | nstart too low | Set nstart = 50-100 |
