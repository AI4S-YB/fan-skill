# Visualization: Time-Series Expression Plots

**Goal:** Generate publication-ready visualizations for time-series
expression data: heatmaps, line plots, cluster profiles, and combined panels.

**Best for:** Processed and clustered time-series expression data.

## Prerequisites

- R 4.0+, `ggplot2`, `pheatmap`, `RColorBrewer`, `reshape2`
- Smoothed expression matrix with cluster/tau assignments
- Time-point metadata

## Heatmap + Line Plot (>= 5 time points)

```r
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

# A. Heatmap with cluster annotation
ann_colors <- list(
  cluster = setNames(
    brewer.pal(n = length(unique(cluster_assignments$cluster)),
               name = "Set1"),
    sort(unique(cluster_assignments$cluster))
  )
)

ann_row <- data.frame(
  cluster = factor(cluster_assignments$cluster),
  row.names = cluster_assignments$gene
)

# Sort genes by cluster, then by within-cluster order
gene_order <- cluster_assignments$gene[order(cluster_assignments$cluster)]

pdf("heatmap_timeseries.pdf", width = 10, height = 12)
pheatmap(
  z_mat[gene_order, ],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main = "Time-Series Expression Heatmap",
  gaps_row = cumsum(table(cluster_assignments$cluster)),
  labels_col = time_labels
)
dev.off()

# B. Cluster center line plots
plot_data <- data.frame()
for (i in 1:nrow(centers)) {
  plot_data <- rbind(plot_data, data.frame(
    cluster = factor(i),
    time = time_points,
    expression = centers[i, ],
    stringsAsFactors = FALSE
  ))
}

pdf("cluster_centers_lines.pdf", width = 10, height = 8)
ggplot(plot_data, aes(x = time, y = expression,
                       color = cluster, group = cluster)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Time point", y = "Standardized expression",
       title = "Cluster Center Profiles") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")
dev.off()

# C. Combined: small multiples with cluster centers + individual genes
pdf("cluster_multiples.pdf", width = 14, height = 10)
par(mfrow = c(ceiling(k/3), 3), mar = c(3, 3, 2, 1))
for (i in 1:k) {
  cluster_genes <- names(km$cluster[km$cluster == i])
  cluster_mat <- z_mat[cluster_genes, ]

  # Plot all individual genes in transparent grey
  matplot(time_points, t(cluster_mat),
          type = "l", lty = 1, col = rgb(0, 0, 0, 0.1),
          xlab = "", ylab = "", main = paste("Cluster", i),
          ylim = range(z_mat))
  # Overlay mean in thick color
  lines(time_points, centers[i, ], lwd = 3, col = i + 1)
  abline(h = 0, lty = 2, col = "grey50")
}
dev.off()

# D. Tau specificity dot plot
tau_df <- data.frame(
  gene = names(tau_values),
  tau = tau_values,
  max_tissue = max_tissue[names(tau_values)]
)

pdf("tau_specificity.pdf", width = 12, height = 6)
ggplot(tau_df, aes(x = max_tissue, y = tau, fill = max_tissue)) +
  geom_violin(alpha = 0.5, draw_quantiles = c(0.5)) +
  geom_jitter(alpha = 0.1, width = 0.2, size = 0.5) +
  geom_hline(yintercept = c(0.5, 0.8), lty = 2, color = c("orange", "red")) +
  labs(x = "Tissue / Stage", y = "Tau Index",
       title = "Expression Specificity by Tissue") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()
```

## Simple Line Plot (< 5 time points)

```r
# Simple line plot with mean +/- SE per cluster
library(ggplot2)

# Per-cluster mean and SE for line plot
line_data <- data.frame()
for (i in 1:k) {
  cluster_genes <- names(km$cluster[km$cluster == i])
  cluster_mat <- z_mat[cluster_genes, , drop = FALSE]

  for (tp in 1:ncol(cluster_mat)) {
    line_data <- rbind(line_data, data.frame(
      cluster = factor(i),
      time = time_points[tp],
      mean = mean(cluster_mat[, tp]),
      se = sd(cluster_mat[, tp]) / sqrt(nrow(cluster_mat)),
      n = nrow(cluster_mat),
      stringsAsFactors = FALSE
    ))
  }
}

pdf("simple_line_plot.pdf", width = 8, height = 6)
ggplot(line_data, aes(x = time, y = mean,
                       color = cluster, group = cluster)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se,
                  fill = cluster), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Time point", y = "Standardized expression",
       title = "Time-Series Expression by Cluster",
       subtitle = paste0("Mean +/- SE (", nrow(z_mat), " genes, ", k, " clusters)")) +
  theme_bw(base_size = 14)
dev.off()
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| color palette | RColorBrewer Set1 | Color palette for clusters |
| heatmap scale | row z-score | Standardization for heatmap |
| transparency (alpha) | 0.1 | Alpha for individual gene lines in multiples |
| PDF dimensions | 10x12 | Default PDF size |
| time_labels | column names | X-axis labels for time points |

## Plant-Specific Notes

- **Cluster naming**: label clusters by their trend rather than a number
  (e.g., "early up, late down" instead of "Cluster 3"). Makes figures
  self-explanatory for plant biology audiences.
- **Time-axis labeling**: use hours for circadian, days for developmental
  series, and explicit stages (e.g., "V3", "R1") for field phenology.
  Never use arbitrary index numbers.
- **Multi-panel for tissues**: if you have both a time course and tissue
  comparison, put time-series heatmaps and tau violin plots in the same
  figure panel (2x2 grid). This is the standard in plant transcriptome atlases.
- **Color conventions**: green = up/photosynthetic, red/brown = stress/senescence,
  blue = circadian/clock, consistent with plant journal conventions.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| PDF blank/empty | `dev.off()` not called | Ensure paired `pdf()` / `dev.off()` |
| Heatmap too dense to read | Too many genes plotted | Subsample or aggregate |
| Colors indistinguishable | Too many clusters for Set1 | Use `viridis` palette for > 9 clusters |
| Line plot overly noisy | Too many individual gene traces | Use mean +/- SE ribbon only |
| `ggplot: object not found` | Variable name mismatch in aes() | Check spelling; use `aes_string()` for dynamic variables |
