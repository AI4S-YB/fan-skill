# Typography for Scientific Figures

## Font Family

### Primary choice: Sans-serif

```r
# Default sans-serif (works on all platforms without additional fonts)
base_family = "sans"

# Platform-specific fallbacks:
# - Windows: "Arial"
# - macOS: "Helvetica"
# - Linux: "Liberation Sans" or "Nimbus Sans"
```

Sans-serif fonts are standard for scientific figures because:
- Higher legibility at small sizes than serif fonts
- Cleaner appearance with data ink (less visual clutter)
- Consistent rendering across PDF/SVG/PNG output

### When serif is acceptable

Only use serif fonts (Times New Roman, Computer Modern) when:
- The journal explicitly requires them
- The figure contains substantial text annotations that read like body text
- The figure is embedded in a LaTeX document using Computer Modern

```r
# Computer Modern (requires 'extrafont' or 'showtext' package)
library(showtext)
font_add("cm", "cmunrm.ttf")  # Computer Modern
showtext_auto()
base_family = "cm"
```

## Font Size Hierarchy

All sizes in points (pt). Based on `theme_plant_scientific(base_size = 8)`:

| Element | Size | Justification |
|---------|------|---------------|
| Axis text (tick labels) | 8 pt (base_size) | Most abundant text; must be readable but not dominant |
| Axis title | 10 pt (base_size + 2) | One per axis; needs distinction from tick labels |
| Facet strip text | 9 pt (base_size + 1), bold | Navigation labels for multi-panel; bold for hierarchy |
| Legend text | 8 pt (base_size) | Less critical; should not compete with data |
| Legend title | 9 pt (base_size + 1) | Slightly larger for hierarchy |
| Panel tag (A, B, C...) | 12 pt, bold | Entry point for readers scanning the figure |
| Annotation text | 6-7 pt | Smallest allowed; for supplementary labels like gene names |
| In-plot statistical annotations | 6-7 pt | p-values, significance stars; less critical |

### Minimum size rule

**Never go below 6 pt.** At 6 pt, text is at the edge of legibility when printed. Most journals will reject figures with text below 6 pt.

### Adjusting for figure size

For speakers/presentations (slides): increase all sizes by 2-4 pt.
For multi-panel figures: maintain the hierarchy but be consistent across panels.

```r
# Single-column figure (typical ~3.5 inches wide)
theme_plant_scientific(base_size = 8)

# Double-column / full-page figure (~7 inches wide)
theme_plant_scientific(base_size = 9)

# Presentation slide (16:9, viewed from distance)
theme_plant_scientific(base_size = 14)
```

## Spacing Guidelines

```r
theme(
  # Margin between plot area and panel border
  plot.margin = margin(10, 10, 10, 10, unit = "pt"),

  # Spacing between faceted panels
  panel.spacing = unit(0.8, "lines"),   # ~9.6 pt at base 12pt line height

  # Tick mark length
  axis.ticks.length = unit(0.05, "cm"), # ~1.4 pt — subtle

  # Legend key size
  legend.key.size = unit(0.4, "cm"),    # ~11 pt — large enough to see color

  # Space between legend items
  legend.spacing.y = unit(0.1, "cm"),

  # Margin around legend
  legend.margin = margin(0, 0, 0, 0, unit = "pt")
)
```

## Text Angle Guidelines

| Context | Angle | Justification |
|---------|-------|---------------|
| X-axis labels (2-5 categories) | 0 (horizontal) | Always preferred |
| X-axis labels (6-15 categories) | 45 degrees | Still readable; aligns with tick marks |
| X-axis labels (15+ categories) | 90 degrees | Only when necessary; consider horizontal bar chart instead |
| Y-axis labels | 90 (vertical) | Standard |
| Annotation in narrow space | 90 | Legible if short (<6 characters) |

```r
# 45-degree labels with proper alignment
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# 90-degree labels
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

## Mathematical Notation in ggplot2

Use `expression()` and `bquote()` for proper mathematical typesetting:

```r
# Subscripts and superscripts
labs(x = expression(log[10](P-value)))
labs(y = expression(italic(K)[s]))     # italic K with subscript s
labs(y = expression(-log[10](italic(P))))  # -log10(P) in italics

# Greek letters
labs(x = expression(alpha))
labs(x = expression(Delta*Delta*Ct))   # Delta Delta Ct

# Combining with text and units
labs(y = expression("Plant Height" ~ (cm) * " " * bar(x) %+-% SEM))

# Using bquote for dynamic values
r2 <- 0.85
labs(subtitle = bquote(R^2 == .(round(r2, 2))))
```

## Using Custom Fonts with showtext

```r
library(showtext)

# Add Google fonts
font_add_google("Roboto", "roboto")
font_add_google("Lato", "lato")
font_add_google("Fira Sans", "fira")

# Or add system fonts
font_add("Helvetica", "/System/Library/Fonts/Helvetica.ttc")

# Enable showtext rendering
showtext_auto()

# Use in ggplot2
theme_plant_scientific(base_family = "roboto")
```

## Font Rendering Tips

1. **PDF output**: Use `device = cairo_pdf` in `ggsave()` for proper font embedding. Default `pdf()` may substitute fonts.
2. **PNG output**: Use `device = "ragg"` or `type = "cairo"` in `ggsave()` for anti-aliased text.
3. **Font embedding in PDFs**: Use `extrafont::embed_fonts()` to ensure fonts are embedded for journal submission.

```r
# Best practice for saving with proper font rendering
ggsave("figure.pdf", plot = p, device = cairo_pdf, width = 8, height = 6)
ggsave("figure.png", plot = p, type = "cairo", width = 8, height = 6, dpi = 300)
```
