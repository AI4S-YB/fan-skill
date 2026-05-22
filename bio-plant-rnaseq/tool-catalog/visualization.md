# RNA-seq Visualization

**Goal:** Volcano plot, MA plot, heatmap, PCA, enrichment dot plot
**Best for:** All RNA-seq results

## Volcano Plot

```r
library(ggplot2)
library(ggrepel)

res$sig <- "NS"
res$sig[res$padj < 0.05 & res$log2FoldChange > 1] <- "UP"
res$sig[res$padj < 0.05 & res$log2FoldChange < -1] <- "DOWN"

top_genes <- rownames(res[order(res$padj), ])[1:10]

ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 0.5) +
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red", "NS" = "grey")) +
  geom_text_repel(data = res[top_genes, ], aes(label = top_genes), size = 3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot")
```

## DEG Count Bar Plot

```r
deg_counts <- data.frame(
  contrast = names(deg_lists),
  UP = sapply(deg_lists, function(x) sum(x$log2FoldChange > 1 & x$padj < 0.05)),
  DOWN = sapply(deg_lists, function(x) sum(x$log2FoldChange < -1 & x$padj < 0.05))
)
```

## WGCNA Module Heatmap

Use `pheatmap` with module assignments and trait correlations.

## Plant-Specific Figure Notes

- Multi-panel figure: Volcano + heatmap of top DEGs + enrichment dot plot
- For polyploids, optionally label homeolog pairs in volcano plot
