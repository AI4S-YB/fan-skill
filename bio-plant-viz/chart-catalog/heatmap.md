# Heatmap

**For**: Gene expression matrices, methylation levels, metabolite abundance, trait correlations

## Data Format Required

| Column | Type | Description |
|--------|------|-------------|
| gene_id | character | Row identifier (gene, metabolite, etc.) |
| sample_1...sample_n | numeric | Values to visualize (expression, abundance, etc.) |

Additional metadata files for annotations:
- Column annotation: sample groups (treatment, tissue, timepoint, etc.)
- Row annotation: gene clusters, functional categories, etc.

## R Code (pheatmap)

```r
library(pheatmap)
library(RColorBrewer)
source("theme/theme_plant_scientific.R")

# Load expression matrix
# expr_mat <- read.table("expression_matrix.txt", header = TRUE, row.names = 1)

# Scale by row (z-score for genes across samples)
# expr_scaled <- t(scale(t(expr_mat)))

# Column annotations
# col_anno <- data.frame(
#   Treatment = c("Control", "Control", "Control", "Drought", "Drought", "Drought"),
#   Tissue = c("Root", "Root", "Shoot", "Root", "Root", "Shoot"),
#   row.names = colnames(expr_mat)
# )

# Annotation colors
# anno_colors <- list(
#   Treatment = c(Control = okabe_ito[2], Drought = okabe_ito[6]),
#   Tissue = c(Root = okabe_ito[3], Shoot = okabe_ito[4])
# )

# ---- Heatmap ----
pheatmap(
  expr_scaled,
  scale = "none",               # already z-scored
  clustering_method = "ward.D2",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  annotation_col = col_anno,
  annotation_colors = anno_colors,
  annotation_names_col = TRUE,
  show_rownames = (nrow(expr_scaled) <= 100),
  fontsize_row = if(nrow(expr_scaled) <= 50) 6 else 4,
  fontsize_col = 8,
  fontsize = 8,
  border_color = NA,
  treeheight_row = 15,
  treeheight_col = 15,
  cutree_rows = 3,              # cut row dendrogram into k clusters
  main = "",
  filename = "heatmap.pdf",
  width = 8,
  height = 10
)
```

## R Code (ComplexHeatmap — for complex annotations)

```r
library(ComplexHeatmap)
library(circlize)

# Color mapping
heatmap_colors <- colorRamp2(
  breaks = c(-2, 0, 2),
  colors = c(okabe_ito[2], "white", okabe_ito[6])
)

# Column annotation
col_ha <- HeatmapAnnotation(
  Treatment = col_anno$Treatment,
  Tissue = col_anno$Tissue,
  col = list(
    Treatment = c("Control" = okabe_ito[2], "Drought" = okabe_ito[6]),
    Tissue = c("Root" = okabe_ito[3], "Shoot" = okabe_ito[4])
  ),
  annotation_name_side = "left"
)

# Row annotation (optional, e.g., gene clusters)
# row_ha <- rowAnnotation(Cluster = anno_block(gp = gpar(fill = okabe_ito[1:3])))

# Main heatmap
Heatmap(
  expr_scaled,
  name = "Z-score",
  col = heatmap_colors,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = "pearson",
  show_row_names = (nrow(expr_scaled) <= 50),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = col_ha,
  row_split = 3,
  border = FALSE,
  row_title = NULL,
  heatmap_legend_param = list(direction = "vertical")
)
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `show_rownames` | TRUE if n <= 100 | Too many gene names make heatmap illegible |
| `fontsize_row` | 6 for <=50, 4 for more | Reduce as gene count increases |
| `cutree_rows` | 3 | Number of gene clusters to highlight |
| Scaling | z-score by row | Essential for comparing genes with different baseline expression |
| Distance metric | correlation | Better for expression patterns than Euclidean |
| Color palette | RdBu | Red-white-blue; use viridis for continuous non-diverging data |

## Plant-Specific Notes

- For organ/tissue atlases: order columns by developmental stage, not alphabetically
- For stress time series: order columns chronologically; consider loess-smoothed line plots as an alternative to heatmap
- For multi-species comparison: use different species as column annotation, not as different rows
- Gene IDs: use species-specific naming (Os for rice, At for Arabidopsis, Zm for maize)
- If >500 genes: consider showing only top variable genes; full heatmap goes to supplementary
- For WGCNA module eigengenes: use heatmap as summary, not the full expression matrix
