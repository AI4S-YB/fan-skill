# Bar Chart, Scatter Plot, and Line Chart

**For**: General-purpose visualizations — bar charts for counts/comparisons, scatter plots for correlations, line charts for time series or dose-response

## Part 1: Bar Chart

### Data Format Required

| Column | Type | Description |
|--------|------|-------------|
| Category | character/factor | Group or condition |
| Value | numeric | Height of bar |
| SE | numeric | Standard error or SD |

### ggplot2 Code

```r
library(ggplot2)
library(dplyr)
source("theme/theme_plant_scientific.R")

# Load data
# bar_data <- read.table("bar_data.txt", header = TRUE)

# ---- Bar Chart with Error Bars ----
ggplot(bar_data, aes(x = Category, y = Value, fill = Category)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = Value - SE, ymax = Value + SE),
    width = 0.2,
    size = 0.4
  ) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  labs(
    x = "",
    y = "Value (units)",
    caption = "Error bars: SE"
  ) +
  theme_plant_scientific()

# ---- Grouped Bar Chart ----
ggplot(bar_data, aes(x = Category, y = Value, fill = Subgroup)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = Value - SE, ymax = Value + SE),
    position = position_dodge(width = 0.8),
    width = 0.2,
    size = 0.4
  ) +
  scale_fill_manual(values = okabe_ito) +
  labs(
    x = "",
    y = "Value (units)",
    fill = "Subgroup",
    caption = "Error bars: SE"
  ) +
  theme_plant_scientific()
```

## Part 2: Scatter Plot

### Data Format Required

| Column | Type | Description |
|--------|------|-------------|
| X | numeric | Independent variable |
| Y | numeric | Dependent variable |
| Group | character/factor | Color/symbol grouping |
| Label | character | Labels for key points (optional) |

### ggplot2 Code

```r
# ---- Basic Scatter Plot ----
ggplot(scatter_data, aes(x = X, y = Y)) +
  geom_point(size = 2, alpha = 0.7, color = okabe_ito[3]) +
  labs(
    x = "X Variable (units)",
    y = "Y Variable (units)"
  ) +
  theme_plant_scientific()

# ---- Scatter with Regression Line ----
ggplot(scatter_data, aes(x = X, y = Y)) +
  geom_point(aes(color = Group), size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", size = 0.6) +
  scale_color_manual(values = okabe_ito) +
  labs(
    x = "X Variable (units)",
    y = "Y Variable (units)",
    color = "Group"
  ) +
  theme_plant_scientific()

# ---- Scatter with LOESS Smoothing ----
ggplot(scatter_data, aes(x = X, y = Y)) +
  geom_point(size = 1.5, alpha = 0.5, color = "grey60") +
  geom_smooth(method = "loess", se = TRUE, color = okabe_ito[6], size = 0.8, span = 0.75) +
  labs(
    x = "Time / Dose (units)",
    y = "Response (units)"
  ) +
  theme_plant_scientific()

# ---- Correlation Scatter with Marginal Distributions ----
library(ggExtra)
p <- ggplot(scatter_data, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = okabe_ito) +
  theme_plant_scientific()

ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
```

## Part 3: Line Chart (Time Series / Dose Response)

### Data Format Required

Long format:

| Column | Type | Description |
|--------|------|-------------|
| Time | numeric | Time point, dose, or ordinal position |
| Value | numeric | Measurement value |
| Group | character/factor | Line grouping |
| SE | numeric | Standard error at each point |

### ggplot2 Code

```r
# ---- Line Chart with Error Ribbons ----
ggplot(line_data, aes(x = Time, y = Value, color = Group, fill = Group)) +
  geom_line(size = 0.8) +
  geom_ribbon(
    aes(ymin = Value - SE, ymax = Value + SE),
    alpha = 0.15,
    color = NA
  ) +
  geom_point(size = 2) +
  scale_color_manual(values = okabe_ito) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  labs(
    x = "Time (days)",
    y = "Response (units)",
    color = "Treatment"
  ) +
  theme_plant_scientific()

# ---- Line Chart with Facets for Multiple Traits ----
ggplot(line_data, aes(x = Time, y = Value, color = Group, fill = Group)) +
  geom_line(size = 0.8) +
  geom_ribbon(
    aes(ymin = Value - SE, ymax = Value + SE),
    alpha = 0.15,
    color = NA
  ) +
  geom_point(size = 1.5) +
  scale_color_manual(values = okabe_ito) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  facet_wrap(~ Trait, scales = "free_y", ncol = 2) +
  labs(
    x = "Time (days)",
    y = "Value",
    color = "Treatment"
  ) +
  theme_plant_scientific()

# ---- Individual Trajectories (for repeated measures, show each subject) ----
ggplot(line_data, aes(x = Time, y = Value, group = Subject)) +
  geom_line(alpha = 0.3, size = 0.3, color = "grey50") +
  stat_summary(aes(group = Group, color = Group), fun = mean, geom = "line", size = 1.2) +
  stat_summary(
    aes(group = Group, color = Group),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.5,
    size = 0.5
  ) +
  scale_color_manual(values = okabe_ito) +
  labs(
    x = "Time (days)",
    y = "Response (units)",
    color = "Treatment"
  ) +
  theme_plant_scientific()
```

## Key Parameters

### Bar Chart
| Parameter | Default | Guidance |
|-----------|---------|----------|
| Bar width | 0.7 | Reduce for more groups |
| `position_dodge` width | 0.8 | Slightly wider than bar width for spacing |

### Scatter Plot
| Parameter | Default | Guidance |
|-----------|---------|----------|
| `geom_smooth` method | lm/loess | lm for linear; loess for nonlinear; gam for complex |
| Point transparency | 0.7 | Lower for >500 points |
| `span` (loess) | 0.75 | Lower for more local flexibility; higher for smoother |

### Line Chart
| Parameter | Default | Guidance |
|-----------|---------|----------|
| `geom_ribbon` alpha | 0.15 | Deliberately transparent to show overlap |
| Error bar type | SE or 95% CI | More conservative for small n |
| Point markers | YES | Add points at each time/dose for clarity; remove if very dense |

## Plant-Specific Notes

- Growth curves: use time-series line chart. Plant growth is sigmoid; LOESS or nonlinear model fits better than linear
- Diurnal cycles: time on x-axis in 24h format. Add shaded regions for night periods with `annotate("rect")`
- Dose-response: typically log10 dose on x-axis; fit 4-parameter log-logistic model (drc package)
- Photosynthesis light-response curves: PPFD on x-axis, A (net assimilation) on y-axis; fit rectangular hyperbola
- Multi-environment trials: combine bar chart (for means) and scatter (for individual environment values) in a two-panel figure
- Phenotypic correlation: scatter matrix (pairs or GGally::ggpairs) for >3 traits; do not use separate scattered panels
