# Multi-Panel Figure Composition

**For**: Combining multiple ggplot2 plots into one publication-quality figure using patchwork or cowplot

## When to Combine Plots

A multi-panel figure tells one complete story. Each panel is a component of that story:
- Panel A: Overview / home view (e.g., full Manhattan plot, PCA, or phenotype distribution)
- Panel B: Detail / zoom view (e.g., regional Manhattan, candidate gene expression, or LD block)
- Panel C: Validation / support view (e.g., QQ plot, enrichment analysis, or haplotype analysis)

## patchwork (Recommended)

```r
library(ggplot2)
library(patchwork)
source("theme/theme_plant_scientific.R")

# ---- Basic Composition ----
# Horizontal: A | B
combined <- plot_a + plot_b

# Vertical: A / B
combined <- plot_a / plot_b

# Grid: A | B / C for 2-column layout (C spans both columns)
combined <- (plot_a | plot_b) / plot_c + plot_layout(heights = c(1, 0.8))

# ---- Complex Layout with Annotation ----
combined <- (plot_a | plot_b) / (plot_c | plot_d) +
  plot_annotation(
    tag_levels = "A",          # Auto-label A, B, C, D
    tag_suffix = "",           # No suffix (A vs A))
    theme = theme_plant_scientific()
  ) &
  theme(
    plot.tag = element_text(size = 12, face = "bold")
  )

# ---- Fine Control with plot_layout ----
combined <- (plot_a + plot_b) / plot_c +
  plot_layout(
    ncol = 2,
    nrow = 2,
    widths = c(1, 1),          # equal width columns
    heights = c(1, 0.6),       # plot_c shorter than a+b
    guides = "collect"         # collect duplicate legends
  ) +
  plot_annotation(
    tag_levels = "A",
    caption = "Figure caption text here"
  )

# ---- Inset (small plot inside large plot) ----
main_plot + inset_element(
  small_plot,
  left = 0.6, bottom = 0.55, right = 0.95, top = 0.95
)

# ---- Save Multi-Panel Figure ----
ggsave(
  "figure_combined.pdf",
  plot = combined,
  width = 8,
  height = 10,
  dpi = 600
)

ggsave(
  "figure_combined.png",
  plot = combined,
  width = 8,
  height = 10,
  dpi = 300
)
```

## cowplot (Alternative, Good for Alignment)

```r
library(cowplot)

# ---- Grid Layout ----
combined <- plot_grid(
  plot_a, plot_b, plot_c, plot_d,
  labels = c("A", "B", "C", "D"),
  label_size = 10,
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1),
  rel_heights = c(1, 0.8),
  align = "hv",
  axis = "lr"
)

# ---- Shared Legend ----
# Remove legend from individual plots
legend <- get_legend(plot_a + theme(legend.position = "bottom"))
plots_no_legend <- plot_grid(
  plot_a + theme(legend.position = "none"),
  plot_b + theme(legend.position = "none"),
  ncol = 2
)
combined <- plot_grid(plots_no_legend, legend, ncol = 1, rel_heights = c(1, 0.1))

# ---- Save ----
save_plot("figure_combined.pdf", combined, base_width = 8, base_height = 10)
```

## Common Layout Recipes

### GWAS Figure
```
|---- A: Manhattan (full width) ----|
|-- B: QQ (left) --|-- C: Regional (right)--|
|---- D: LD Heatmap (full width) ----|
```

```r
gwas_fig <- (manhattan_plot) /
             (qq_plot | regional_plot) /
             (ld_heatmap) +
  plot_layout(heights = c(1.2, 1, 0.8)) +
  plot_annotation(tag_levels = "A")
```

### RNA-seq DEG Figure
```
|-- A: Volcano --|-- B: Heatmap (top DEGs) --|
|------------ C: GO Enrichment Dot ------------|
```

```r
deg_fig <- (volcano_plot | heatmap_plot) / enrichment_dot +
  plot_layout(heights = c(1, 0.7)) +
  plot_annotation(tag_levels = "A")
```

### Population Structure Figure
```
|-- A: PCA --|-- B: ADMIXTURE bar --|
|-------- C: NJ Tree ---------------|
```

```r
pop_fig <- (pca_plot | admixture_plot) / nj_tree +
  plot_layout(heights = c(1, 0.8))
```

### GS Results Figure
```
|-- A: Pred vs Obs Scatter --|-- B: GEBV Ranking (top 20) --|
|---------- C: Prediction Accuracy by Model ----------------|
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `tag_levels` | "A" | Standard uppercase; use "a" for lowercase; "1" for numbers |
| `guides = "collect"` | TRUE | Use when multiple panels share the same color legend |
| `align` (cowplot) | "hv" | Align both horizontal and vertical |
| `rel_widths` | c(1, 1) | Adjust for panels needing more space (e.g., tree plot wider) |
| `rel_heights` | c(1, 1) | Adjust for panels needing more height (e.g., heatmap taller) |

## Best Practices

1. **One legend for all panels.** Use `guides = "collect"` in patchwork or extract shared legend with cowplot.
2. **Consistent theme across all panels.** Apply `theme_plant_scientific()` to each individual plot; the combined theme inherits.
3. **Panel labels (A, B, C...) always upper-left or upper-right.** Be consistent. patchwork places them top-left by default.
4. **Each panel should be self-explanatory.** A reader looking at any single panel should understand what's being shown.
5. **Don't exceed 6 panels per figure.** More than 6 panels becomes visually overwhelming. Split into multiple figures.
6. **Align axes when panels share dimensions.** Use `align = "hv"` or match ylim/xlim explicitly.

## Plant-Specific Notes

- Figures with plant photographs (phenotype, histology): use cowplot `draw_image()` or `ggdraw()` + `draw_image()` to add images as panels
- Chromosome ideograms: can be combined with GWAS results using custom grobs or karyoploteR output
- For multi-species figures: add species labels as strip annotation or panel titles, not in the legend
- Publication-specific: some plant journals (Plant Cell, Plant Physiology, New Phytologist) have specific figure width requirements; check before finalizing layout
