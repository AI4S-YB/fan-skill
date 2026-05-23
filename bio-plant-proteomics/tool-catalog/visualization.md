# Proteomics Visualization

**Goal:** Volcano plots, heatmaps, PCA, PPI network visualizations for plant proteomics
**Best for:** Publication-quality figures and exploratory analysis

## Volcano Plot

```r
library(ggplot2)
library(ggrepel)

results <- read.csv("dep_results.csv")

results$category <- "NS"
results$category[results$adj.P.Val < 0.05 & results$logFC > 1] <- "Up"
results$category[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Down"

# Top proteins to label
top10 <- head(results[order(results$adj.P.Val), ], 10)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = category), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("Up" = "#e74c3c",
                                 "Down" = "#3498db",
                                 "NS" = "grey80")) +
  geom_text_repel(data = top10,
                  aes(label = Gene.names),
                  size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  theme_minimal(base_size = 12) +
  labs(title = "Differential Protein Abundance",
       subtitle = paste0("Up: ", sum(results$category == "Up"),
                         " | Down: ", sum(results$category == "Down")),
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value")
ggsave("volcano_plot.png", width = 8, height = 6, dpi = 300)
```

## Heatmap of DEPs

```r
library(pheatmap)

# Select significantly changed proteins
sig_proteins <- results[results$significant, "Protein.IDs"]
sig_mat <- lfq_imp[rownames(lfq_imp) %in% sig_proteins, ]

# Z-score normalization
sig_mat_z <- t(scale(t(sig_mat)))

# Annotation
annotation_col <- data.frame(
  Condition = groups,
  row.names = colnames(sig_mat_z)
)

pheatmap(sig_mat_z,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         scale = "none",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("#3498db", "white", "#e74c3"))(100),
         main = "Differentially Abundant Proteins",
         fontsize = 10,
         cutree_rows = 3)
```

## PCA Plot

```r
library(ggfortify)

# PCA on all proteins
pca_res <- prcomp(t(lfq_imp), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Condition <- groups

# Variance explained
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  labs(x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"),
       title = "PCA — Protein Abundance") +
  theme_minimal()
```

## Phosphorylation Motif Logo

```python
import logomaker
import pandas as pd

# Extract 15-mer centered on phosphosite
sequences = phos_filtered["window"].tolist()

# Create position weight matrix
from collections import Counter
pwm = pd.DataFrame(0.0, index=list("ACDEFGHIKLMNPQRSTVWY"),
                   columns=range(-7, 8))

for seq in sequences:
    if len(seq) == 15:
        for i, aa in enumerate(seq):
            pos = i - 7
            if aa in pwm.index:
                pwm.loc[aa, pos] += 1

# Normalize by column
pwm_norm = pwm / pwm.sum()

# Create logo
fig, ax = plt.subplots(figsize=(12, 4))
logomaker.Logo(pwm_norm, ax=ax)
ax.set_title("Phosphorylation Site Motif (±7 residues)")
ax.set_xlabel("Position relative to phosphosite")
plt.savefig("phospho_motif_logo.png", dpi=300, bbox_inches='tight')
```

## PPI Network Visualization in Cytoscape

### Styling the Network

1. Import STRING network into Cytoscape
2. Map node color: `logFC` (red = up, blue = down, white = NS)
3. Map node size: `-log10(adj.P.Val)` (larger = more significant)
4. Map edge thickness: `combined_score` (thicker = higher confidence)
5. Apply layout: `Prefuse Force Directed` or `yFiles Organic`

### Export

```
File > Export > Network to Image
Format: PDF or PNG
Resolution: 300 DPI
Export text as font: Yes (for PDF)
```

## Multi-Panel Figure Assembly

Combine plots using `patchwork` in R:

```r
library(patchwork)

(volcano | heatmap) / (pca | bar) +
  plot_annotation(
    title = "Plant Proteomics Analysis",
    tag_levels = 'a'
  )
```

## Plant-Specific Figure Notes

- Always label RuBisCO large/small subunits on volcano plots (they dominate the plot)
- For phosphoproteomics: use separate color scheme (orange = phospho-up, purple = phospho-down)
- PPI networks often cluster by subcellular location: plastid, mitochondria, cytosol, nucleus
- Time-series phosphoproteomics: use line plots showing temporal phosphorylation dynamics
- Include QC metrics in supplementary: peptide count distribution, mass error distribution, missed cleavages

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Volcano "shark fin" shape | Non-logged or poorly normalized data | Log2 transform and normalize |
| PCA no separation | Genuine lack of difference or batch effect dominates | Check batch covariates, include in model |
| Heatmap all one color | Zero-variance genes dominate after filtering | Remove low-variance proteins before clustering |
| PPI network "hairball" | Too many interactions with low confidence | Filter to high confidence (>= 700) |
