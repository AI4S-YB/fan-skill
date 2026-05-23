# Metabolomics Visualization

**Goal:** Generate publication-quality figures for plant metabolomics: PCA, volcano plots, heatmaps, pathway views, and compound class distributions
**Best for:** All metabolomics studies — essential for data exploration and publication

## Prerequisites
- R 4.0+ with ggplot2, pheatmap, plotly, pathview
- Python 3 with matplotlib, seaborn (alternative)

## PCA Scores Plot

```r
library(ggplot2)

# PCA on log2-transformed feature table
pca <- prcomp(t(feature_table_log2), center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca$x)
pca_scores$group <- metadata$group
pca_scores$sample <- rownames(pca_scores)

# Variance explained
var_exp <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = group, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  labs(
    title = "PCA Scores Plot",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Set1")
```

## Volcano Plot

```r
library(ggplot2)
library(ggrepel)

results <- read.csv("differential_metabolites.csv")
results$neg_log10_p <- -log10(results$P.Value)
results$significant <- results$adj.P.Val < 0.05

# Label top 10 features
top_hits <- results[order(results$P.Value), ][1:10, ]

ggplot(results, aes(x = logFC, y = neg_log10_p, color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_hits,
                  aes(label = putative_name),
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("grey60", "#C73E1D"),
                     labels = c("NS", "Significant (q<0.05)")) +
  labs(x = "log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot: Treatment vs Control") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
```

## Compound Class Distribution (CANOPUS)

```r
library(ggplot2)

# Read CANOPUS results
canopus <- read.csv("canopus_summary.csv")

# Count compounds per class (top level)
class_counts <- as.data.frame(table(canopus$superclass))
names(class_counts) <- c("Class", "Count")
class_counts <- class_counts[order(-class_counts$Count), ]

# Select top 15 classes
top_classes <- head(class_counts, 15)

ggplot(top_classes, aes(x = reorder(Class, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
  coord_flip() +
  labs(x = "", y = "Number of Features",
       title = "Compound Class Distribution (CANOPUS)") +
  theme_bw(base_size = 12)
```

## Heatmap with Hierarchical Clustering

```r
library(pheatmap)
library(RColorBrewer)

# Extract significant features
sig_features <- results[results$adj.P.Val < 0.05, "feature_id"]
sig_matrix <- feature_table_log2[sig_features, ]

# Z-score normalization
sig_matrix_z <- t(scale(t(sig_matrix)))

# Column annotation
annotation_col <- data.frame(
  Group = metadata$group,
  row.names = rownames(metadata)
)

# Color palette
ann_colors <- list(
  Group = c(Control = "#2E86AB", Treatment = "#C73E1D")
)

pheatmap(sig_matrix_z,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
         main = "Significant Metabolites (Z-score, q < 0.05)",
         fontsize = 10,
         border_color = NA)
```

## Boxplot of Top Metabolites

```r
# Plot individual metabolite intensities
top_metabolites <- head(sig_results$feature_id, 9)
plot_data <- data.frame()
for (feat in top_metabolites) {
  temp <- data.frame(
    intensity = as.numeric(feature_table_log2[feat, ]),
    group = metadata$group,
    metabolite = feat
  )
  plot_data <- rbind(plot_data, temp)
}

ggplot(plot_data, aes(x = group, y = intensity, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 3) +
  labs(y = "Log2 Intensity", x = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8))
```

## Interactive Scatter Plot (plotly)

```r
library(plotly)

p <- plot_ly(
  results,
  x = ~logFC,
  y = ~neg_log10_p,
  type = "scatter",
  mode = "markers",
  color = ~significant,
  colors = c("grey60", "#C73E1D"),
  text = ~paste("Feature:", feature_id,
                "<br>m/z:", round(mz, 4),
                "<br>RT:", round(rt, 2), "s",
                "<br>log2FC:", round(logFC, 3),
                "<br>p:", signif(P.Value, 3)),
  hoverinfo = "text",
  marker = list(size = 6, opacity = 0.6)
)
p <- p %>% layout(
  title = "Interactive Volcano Plot",
  xaxis = list(title = "log2 Fold Change"),
  yaxis = list(title = "-log10(p-value)")
)
htmlwidgets::saveWidget(p, "volcano_interactive.html")
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| stat_ellipse level | Confidence ellipse (0.95 for 95% CI) |
| clustering_method | Linkage method for hierarchical clustering |
| colorRampPalette | Color gradient for heatmaps |

## Plant-Specific Visualization Notes

- Compound class distribution pie/bar charts are informative for showing secondary metabolite diversity
- Color-code by pathway: flavonoid (orange), terpenoid (green), alkaloid (purple), phenolic (blue)
- For pathway maps, use plant-specific KEGG codes (ath, osa, zma, sly)
- Side-by-side boxplots for different plant tissues or developmental stages show tissue-specific metabolism
- Interactive plots help explore large metabolomics datasets (1000+ features)

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| PCA not separating groups | Metabolome variation driven by other factors (tissue, batch) | Color by other metadata variables to find confounders |
| Volcano plot asymmetric | Systematic ion suppression in one group | Check extraction and sample preparation |
| Heatmap uninformative | Too many features plotted | Filter to significant only (q < 0.05) |
| Boxplot y-axis misleading | Different baseline intensities | Use Z-score for comparison plots |
