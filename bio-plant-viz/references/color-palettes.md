# Color Palettes for Plant Biology Figures

## Principles

1. **Colorblind-safe by default.** ~8% of males and ~0.5% of females have some form of color vision deficiency (CVD). Use palettes that remain distinguishable under protanopia (red-blind) and deuteranopia (green-blind).
2. **Perceptually uniform.** Steps in the color scale should correspond to equal steps in the data (viridis family).
3. **Grayscale-compatible.** Figures should remain interpretable when printed in black and white. Use luminance differences or add redundant encodings (shape, line type).
4. **Context-appropriate.** Plant biology has natural color associations (green = plant, brown = soil, blue = water, red/yellow = stress).

## Palette 1: Okabe-Ito (Colorblind-Safe, Categorical)

The gold standard for up to 8 categorical groups. Validated for all common forms of colorblindness.

```r
# Okabe-Ito palette (8 colors)
okabe_ito <- c(
  "#E69F00",  # Orange
  "#56B4E9",  # Sky Blue
  "#009E73",  # Bluish Green
  "#F0E442",  # Yellow
  "#0072B2",  # Blue
  "#D55E00",  # Vermillion
  "#CC79A7",  # Reddish Purple
  "#000000"   # Black
)

# Usage in ggplot2
ggplot(data, aes(x = x, y = y, color = group)) +
  scale_color_manual(values = okabe_ito)
```

**Best for**: Group comparisons with <=8 groups (e.g., treatments, genotypes, tissues).

**Avoid**: More than 8 groups -- switch to viridis or use shapes/line types as redundant encodings.

## Palette 2: Viridis (Perceptually Uniform, Sequential)

Four variants, all colorblind-safe and perceptually uniform.

```r
library(ggplot2)

# Viridis D (default, perceptually uniform)
scale_color_viridis_c(option = "D")  # continuous
scale_color_viridis_d(option = "D")  # discrete

# Viridis variants:
# option = "A" : magma   (black-purple-red-yellow)
# option = "B" : inferno (black-dark red-yellow)
# option = "C" : plasma  (dark blue-red-yellow)
# option = "D" : viridis (dark blue-green-yellow)
# option = "E" : cividis (blue-yellow, optimized for CVD)
# option = "F" : rocket  (dark blue-red)
# option = "G" : mako    (dark green-blue)
# option = "H" : turbo   (rainbow — NOT colorblind safe, use only when necessary)

# Reversed direction (darker = more extreme)
scale_fill_viridis_c(option = "D", direction = -1)
```

**Best for**: Heatmaps, continuous color scales, enrichment dot plots.

**When to use which**:
- `viridis (D)`: General purpose. Best overall.
- `inferno (B)`: High dynamic range, emphasizes extremes.
- `cividis (E)`: Maximum colorblind safety.
- `magma (A)`: Bright extremes, dark middle.

## Palette 3: Plant-Specific Palettes

### Field/Soil Palette (Categorical, up to 6 groups)

```r
plant_field <- c(
  "#228B22",  # ForestGreen — healthy plant / control
  "#8B4513",  # SaddleBrown — soil / root
  "#6B8E23",  # OliveDrab  — mature/stressed plant
  "#DAA520",  # Goldenrod  — harvest / yield / grain
  "#556B2F",  # DarkOliveGreen — canopy / leaf
  "#CD853F"   # Peru — dry soil / drought condition
)
```

**Use for**: Field trials, soil treatments, biomass/yield comparisons, canopy measurements.

### Stress Response Palette (Diverging, blue-green-red)

```r
plant_stress <- c(
  "#1B9E77",  # Teal green — control / unstressed
  "#66A61E",  # Light green — mild stress
  "#E6AB02",  # Gold — moderate stress
  "#D95F02",  # Orange — severe stress
  "#A6761D",  # Brown — extreme stress
  "#7570B3"   # Purple — recovery
)
```

**Use for**: Drought, heat, salt, nutrient stress gradients. Diverging from control (green) through stressed (orange/red).

### Temperature Palette (Sequential, cool to hot)

```r
plant_temp <- c(
  "#4575B4",  # Blue — cold
  "#91BFDB",  # Light blue — cool
  "#E0F3F8",  # Very light blue — mild
  "#FEE090",  # Light orange — warm
  "#FC8D59",  # Orange — hot
  "#D73027"   # Red — heat stress
)
```

**Use for**: Temperature response curves, climate gradients, growing degree days.

### Photosynthesis Palette (Sequential, green)

```r
# Light/chlorophyll gradient
plant_photo <- c(
  "#00441B", "#006D2C", "#238B45", "#41AB5D",
  "#74C476", "#A1D99B", "#C7E9C0", "#F7FCF5"
)
```

**Use for**: Chlorophyll content, photosynthetic rate, NDVI, biomass.

## Color Brewer Palettes

RColorBrewer provides additional palettes:

```r
library(RColorBrewer)

# Sequential (9 levels max)
brewer.pal(9, "Greens")     # Plant biomass, chlorophyll
brewer.pal(9, "YlOrBr")     # Drought stress
brewer.pal(9, "PuBu")       # Water potential, precipitation
brewer.pal(9, "YlGn")       # Yield, growth

# Diverging (11 levels max)
brewer.pal(11, "BrBG")      # Brown-green (soil-plant divergence)
brewer.pal(11, "RdBu")      # Heatmap (up-regulated/down-regulated)
brewer.pal(11, "PiYG")      # Purple-green (GO terms, enrichment up/down)
brewer.pal(11, "PuOr")      # Purple-orange

# Qualitative (8-12 levels)
brewer.pal(8, "Set1")       # General purpose, similar to Okabe-Ito
brewer.pal(8, "Dark2")      # Muted — good for secondary groups
brewer.pal(12, "Paired")    # 12 colors, paired design (before/after, homo/heterozygote)
```

## Practical Guidelines

### How to choose a palette

```
Decision tree:

What kind of data?
├── Categorical (groups)
│   ├── 2-8 groups → Okabe-Ito
│   ├── 9-12 groups → RColorBrewer Paired or Set3
│   └── >12 groups → Facet instead, or accept that colors will blend
│
├── Sequential (low to high)
│   ├── General purpose → viridis D
│   ├── Plant-specific context → plant_photo or plant_temp
│   └── Diverging (below/above baseline) → BrBG or RdBu
│
└── Diverging (two directions from center)
    ├── Up/down regulation → RdBu (red = up, blue = down) or viridis
    └── Control vs treatment → BrBG (brown = soil/control, green = plant/treatment)
```

### Testing for colorblindness

```r
# Install colorblindr to simulate CVD in your plots
# remotes::install_github("clauswilke/colorblindr")
library(colorblindr)
cvd_grid(my_plot)  # Shows plot under 4 types of colorblindness
```

### Redundant encodings

When color alone is insufficient (or if grayscale printing is required):

- Add point shapes: `aes(shape = group)` + `scale_shape_manual(values = c(16, 17, 15, 18))`
- Add line types: `aes(linetype = group)` + `scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash"))`
- Add direct labels: `geom_text_repel()` for key data points
- Use facets: `facet_wrap(~ group)` for small multiples
