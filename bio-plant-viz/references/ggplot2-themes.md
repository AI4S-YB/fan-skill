# ggplot2 Theme Reference

## theme_plant_scientific()

### Philosophy

`theme_plant_scientific()` is built on three principles:

1. **Maximize data-ink ratio (Tufte).** Remove non-data elements: top/right borders, minor gridlines, dark backgrounds.
2. **Prioritize readability.** 8pt base font, clear hierarchy, adequate spacing.
3. **Publication-ready by default.** Sensible defaults that work for most plant biology journals. Customize when needed.

### Source Location

The theme is defined in `theme/theme_plant_scientific.R`. Load it with:

```r
source("theme/theme_plant_scientific.R")
```

### Complete Function

```r
theme_plant_scientific <- function(base_size = 8, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Borders: bottom and left only (Tufte)
      panel.border = element_rect(fill = NA, color = "black", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.2, color = "grey90"),

      # Clean axis lines
      axis.line = element_blank(),
      axis.ticks = element_line(size = 0.3),
      axis.ticks.length = unit(0.05, "cm"),

      # Typography hierarchy
      axis.title = element_text(size = base_size + 2, face = "plain"),
      axis.text = element_text(size = base_size),
      strip.text = element_text(size = base_size + 1, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 1),
      legend.key = element_rect(fill = NA, color = NA),
      legend.key.size = unit(0.4, "cm"),

      # Faceting
      panel.spacing = unit(0.8, "lines"),

      # Default aspect ratio
      aspect.ratio = 0.75
    )
}
```

### How It Differs from Built-in Themes

| Feature | theme_bw | theme_classic | theme_minimal | theme_plant_scientific |
|---------|----------|---------------|---------------|----------------------|
| Panel border | Full box | Open bottom+left | None | Full (thin, black) |
| Gridlines | White + grey | None | Grey | Major only (grey90, 0.2pt) |
| Panel background | White | White | White | White |
| Strip background | Grey85 | Grey85 | NA | Grey95, no border |
| Strip text | 0.8*base, plain | 0.8*base, plain | 0.8*base, plain | base+1, bold |
| Legend background | White | NA | NA | NA |
| Base font size | 11 | 11 | 11 | 8 (publication-optimized) |
| Aspect ratio | NA | NA | NA | 0.75 |
| Tick length | default | default | default | 0.05cm (short) |

### Why base_size = 8?

- Most journals print figures at 1-column width (~3.5 inches). At 8pt, ~10 characters fit per cm, making axis labels comfortably readable.
- ggplot2's default 11pt is designed for screen viewing. 8pt is designed for print.
- Increase to 9-10pt for presentations, posters, or 2-column figures.

## Customization Patterns

### Override Specific Elements

```r
# Make a plot for a presentation
my_plot + theme_plant_scientific(base_size = 14)

# Remove gridlines entirely (for clean scatter plots)
my_plot + theme_plant_scientific() +
  theme(panel.grid.major = element_blank())

# Add bottom axis line only (for bar charts)
my_plot + theme_plant_scientific() +
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.3)
  )

# Adjust aspect ratio for specific chart types
# Manhattan plot: very wide
my_manhattan + theme_plant_scientific() + theme(aspect.ratio = 0.3)

# Heatmap: depends on data shape
my_heatmap + theme_plant_scientific() + theme(aspect.ratio = NULL)

# Square scatter (1:1 relationship)
my_scatter + theme_plant_scientific() + theme(aspect.ratio = 1)
```

### Legend Placement

```r
# Top (for 2-4 groups, saves vertical space)
theme(legend.position = "top", legend.direction = "horizontal")

# Inside plot area (for sparse scatter plots, top-left corner)
theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

# Right (default, for many groups)
theme(legend.position = "right")

# Bottom (for multi-panel figures with collected legend)
theme(legend.position = "bottom")

# Remove (when groups are labeled directly)
theme(legend.position = "none")
```

### Axis Text Handling

```r
# Long x-axis labels: rotate 45 degrees
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Scientific notation for y-axis
scale_y_continuous(labels = scales::scientific)

# Percentage format
scale_y_continuous(labels = scales::percent)

# Custom number format
scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))
```

### Multiple Group Palettes

```r
# Use with viridis (continuous + colorblind-safe)
ggplot(data, aes(x, y, color = z)) +
  scale_color_viridis_c() +
  theme_plant_scientific()

# Use with Okabe-Ito (discrete + colorblind-safe)
ggplot(data, aes(x, y, color = group)) +
  scale_color_manual(values = okabe_ito) +
  theme_plant_scientific()
```

## Journal-Specific Adjustments

### Nature family (Nature, Nature Genetics, Nature Plants)
- Single column: 89 mm wide
- Double column: 183 mm wide
- Font: sans-serif, minimum 6pt
- Color: free but must be colorblind-accessible
- Format: vector (PDF/EPS) preferred

```r
theme_plant_scientific(base_size = 7)  # slightly smaller for tight layouts
ggsave("figure.pdf", width = 89, height = 100, units = "mm")
```

### Plant Cell / Plant Physiology
- Single column: ~80 mm wide
- Double column: ~170 mm wide
- Font: Arial or Helvetica, minimum 7pt
- Format: PDF or EPS

```r
theme_plant_scientific(base_size = 8, base_family = "sans")
```

### New Phytologist
- Figure width: 80 mm (single) or 170 mm (double)
- High-contrast figures encouraged
- Color figures free
- Resolution: 300 DPI minimum for halftones; 600 DPI for line art

### Frontiers in Plant Science
- Figures embedded in text, unlimited color
- Recommended width: 85 mm (small), 120 mm (medium), 180 mm (large)
- Resolution: 300 DPI

## Complete Figure Script Template

```r
#!/usr/bin/env Rscript
# Figure X: [Description]
# Generated: [Date]
# Data source: [path]
library(ggplot2)
library(dplyr)
library(patchwork)
source("theme/theme_plant_scientific.R")

# ---- Load data ----
data <- read.table("data.txt", header = TRUE)

# ---- Panel A: [Description] ----
p_a <- ggplot(data, aes(x = x, y = y, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = okabe_ito) +
  labs(x = "X-axis (units)", y = "Y-axis (units)") +
  theme_plant_scientific()

# ---- Panel B: [Description] ----
p_b <- ggplot(data, aes(...)) +
  ... +
  theme_plant_scientific()

# ---- Combine ----
fig <- (p_a | p_b) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

# ---- Export ----
ggsave("figure_X.pdf", plot = fig, device = cairo_pdf, width = 8, height = 5)
ggsave("figure_X.png", plot = fig, width = 8, height = 5, dpi = 300)
ggsave("figure_X.svg", plot = fig, width = 8, height = 5)

# ---- Session Info ----
sessionInfo()
```
