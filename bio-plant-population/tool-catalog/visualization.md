# Population Genetics Visualization

**Goal:** PCA scatter plot, ADMIXTURE bar plot, NJ tree, Fst heatmap
**Best for:** All population genetics results

## PCA Plot

```r
library(ggplot2)
pca <- read.table("pca_result.eigenvec")
eigenval <- scan("pca_result.eigenval")
pct <- round(eigenval/sum(eigenval)*100, 1)

ggplot(pca, aes(x=V3, y=V4)) +
    geom_point(size=2, alpha=0.7) +
    xlab(paste0("PC1 (", pct[1], "%)")) +
    ylab(paste0("PC2 (", pct[2], "%)")) +
    theme_minimal()
```

## ADMIXTURE Bar Plot

Use pophelper R package for publication-quality ADMIXTURE plots.

## Fst Heatmap

```r
library(pheatmap)
fst <- as.matrix(read.table("fst_result.fst", row.names=1))
pheatmap(fst, display_numbers=TRUE, main="Pairwise Fst")
```
