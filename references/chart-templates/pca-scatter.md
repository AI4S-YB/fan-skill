# PCA Scatter Plot

**For**: Population structure, sample clustering, quality control

## Data Format Required

Two options:

### Option A: Raw genotype matrix for PCA computation
| Column | Type | Description |
|--------|------|-------------|
| ID | character | Sample identifier |
| SNP1...SNPn | numeric | Genotype calls (0/1/2 or dosage) |

Note: First column is sample ID; remaining columns are markers. The PCA is computed from this matrix, then visualized.

### Option B: Pre-computed PCA coordinates
| Column | Type | Description |
|--------|------|-------------|
| ID | character | Sample identifier |
| PC1 | numeric | PC1 score |
| PC2 | numeric | PC2 score |
| PC3 | numeric | PC3 score (optional) |
| Population | character/factor | Group/population label |
| Species | character/factor | (optional) for multi-species analysis |

## ggplot2 Code

```r
library(ggplot2)
source("theme/theme_plant_scientific.R")

# ---- Option 1: Compute PCA from genotype matrix ----
# geno: matrix with samples in rows, markers in columns
# pca_result <- prcomp(geno, center = TRUE, scale. = TRUE)
# pca_scores <- as.data.frame(pca_result$x)
# pca_scores$ID <- rownames(pca_scores)

# Variance explained
# var_exp <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# ---- Option 2: Use pre-computed coordinates ----
# pca_scores <- read.table("pca_scores.txt", header = TRUE)
# var_exp <- c(18.5, 12.3)  # from analysis output

# ---- Plot ----
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, size = 0.4, show.legend = FALSE) +
  scale_color_manual(values = okabe_ito) +
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    color = "Population"
  ) +
  theme_plant_scientific() +
  theme(legend.position = "right")

# ---- With sample labels (small datasets < 50 samples) ----
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, size = 0.4, show.legend = FALSE) +
  geom_text_repel(aes(label = ID), size = 2, max.overlaps = 30) +
  scale_color_manual(values = okabe_ito) +
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    color = "Population"
  ) +
  theme_plant_scientific()

# ---- PC1 vs PC3 (for populations poorly separated by PC1-PC2) ----
ggplot(pca_scores, aes(x = PC1, y = PC3, color = Population)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, size = 0.4, show.legend = FALSE) +
  scale_color_manual(values = okabe_ito) +
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC3 (", var_exp[3], "%)"),
    color = "Population"
  ) +
  theme_plant_scientific()

# ---- Multi-panel: PC1-2 + PC1-3 + scree plot (patchwork) ----
# Scree plot
# scree_data <- data.frame(
#   PC = paste0("PC", 1:10),
#   Variance = var_exp[1:10]
# )
# scree_plot <- ggplot(scree_data, aes(x = reorder(PC, -Variance), y = Variance)) +
#   geom_col(fill = okabe_ito[2]) +
#   geom_text(aes(label = paste0(Variance, "%")), vjust = -0.3, size = 2.5) +
#   labs(x = "", y = "Variance Explained (%)") +
#   theme_plant_scientific() +
#   theme(aspect.ratio = 0.5)
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `level` (ellipse) | 0.95 | 95% CI; use 0.68 for 1-SD ellipse |
| Point size | 2.5 | Reduce for >500 samples |
| Point alpha | 0.8 | Reduce for >500 samples to show density |
| Label samples | FALSE | Only for <50 samples or key samples |

## Plant-Specific Notes

- For breeding populations: color by breeding program/subpopulation/pedigree
- For germplasm collections: color by geographic origin or botanical variety
- For polyploids: run PCA on dosage matrix or per-subgenome markers
- If population structure is weak: show first 4 PCs in a pairs plot
- For multi-species analysis: add shape aesthetic (shape = Species) as redundant encoding
- Scree plot helps justify the number of PCs shown — include it as supplementary figure
