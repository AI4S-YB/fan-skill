# Volcano Plot + MA Plot

**For**: Differential expression results (RNA-seq, proteomics)

## Data Format Required

DESeq2/edgeR/limma output:

| Column | Type | Description |
|--------|------|-------------|
| gene_id | character | Gene or transcript identifier |
| log2FC | numeric | log2 fold change (treatment / control) |
| pvalue | numeric | Raw p-value |
| padj | numeric | FDR-adjusted p-value (BH) |
| baseMean | numeric | Mean normalized count (optional, for MA plot) |

## Volcano Plot (ggplot2 + ggrepel)

```r
library(ggplot2)
library(ggrepel)
library(dplyr)
source("theme/theme_plant_scientific.R")

# Load DEG results
# deg <- read.table("deseq2_results.txt", header = TRUE)

# Classify genes
deg <- deg %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FC > 1   ~ "Up",
      padj < 0.05 & log2FC < -1  ~ "Down",
      TRUE                        ~ "NS"
    ),
    neg_log10_padj = -log10(padj)
  )

# Pick top genes to label
top_genes <- deg %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(20)

# Color palette
volcano_colors <- c("Up" = "#D55E00", "Down" = "#0072B2", "NS" = "grey70")

ggplot(deg, aes(x = log2FC, y = neg_log10_padj, color = significance)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = volcano_colors) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", size = 0.3) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_id),
    size = 2.5,
    max.overlaps = 20,
    box.padding = 0.3
  ) +
  labs(
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10](adjusted ~ italic(P) ~ value)),
    color = "DEGs"
  ) +
  theme_plant_scientific() +
  theme(legend.position = c(0.12, 0.85))

# Annotate DEG counts
deg_counts <- deg %>%
  summarise(
    up   = sum(significance == "Up"),
    down = sum(significance == "Down"),
    ns   = sum(significance == "NS")
  )

# Add count annotation as subtitle (passed to labs or separate text grob)
labs(subtitle = paste0("Up: ", deg_counts$up, " | Down: ", deg_counts$down))
```

## MA Plot (ggplot2)

```r
# MA plot: average expression (A) vs log2 fold change (M)
ggplot(deg, aes(x = log10(baseMean), y = log2FC, color = significance)) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = volcano_colors) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey40", size = 0.3) +
  labs(
    x = expression(log[10]("Mean normalized count")),
    y = expression(log[2] ~ "Fold Change"),
    color = "DEGs"
  ) +
  theme_plant_scientific()
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `log2FC` threshold | 1 (2-fold) | Reduce to 0.58 (1.5-fold) for proteomics |
| `padj` threshold | 0.05 | Use 0.01 for more stringent; 0.1 for exploratory |
| Point alpha | 0.6 | Lower for >10000 genes; higher for <1000 |
| `top_genes` n | 20 | Label fewer for cleaner plot; more for detailed view |

## Plant-Specific Notes

- For polyploid species: consider labeling homeologs with prefix (e.g., TraesCS1A_..., TraesCS1B_...)
- If many genes are significant (>2000): reduce top_genes or switch to MA plot which is less crowded
- Plant RNA-seq often has higher fold changes than mammalian; adjust x-axis limits accordingly
- For time-series DEGs: add a faceted volcano plot with `facet_wrap(~ timepoint)` or use smear plot
- Consider showing vst/rlog normalized counts in a separate heatmap panel
