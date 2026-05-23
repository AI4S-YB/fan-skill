# Box Plot + Violin Plot + Jitter Overlay

**For**: Trait distributions, gene expression by group, phenotype comparisons across treatments

## Data Format Required

Long format:

| Column | Type | Description |
|--------|------|-------------|
| Sample | character/factor | Sample identifier |
| Group | character/factor | Treatment, genotype, tissue, etc. |
| Value | numeric | The measurement (trait value, expression, etc.) |
| Color_by | character/factor | Secondary grouping (optional) |

## ggplot2 Code

```r
library(ggplot2)
library(dplyr)
library(ggbeeswarm)  # for geom_quasirandom
source("theme/theme_plant_scientific.R")

# Load data
# data_long <- read.table("trait_data.txt", header = TRUE)

# Compute sample sizes for annotation
sample_n <- data_long %>%
  group_by(Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(Group, " (n=", n, ")"))

# ---- Box Plot ----
ggplot(data_long, aes(x = Group, y = Value)) +
  geom_boxplot(
    aes(fill = Group),
    width = 0.6,
    outlier.shape = NA,    # suppress default outliers (will add jitter)
    alpha = 0.7
  ) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5, color = "grey30") +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  scale_x_discrete(labels = sample_n$label) +
  labs(
    x = "",
    y = "Trait Value (units)"
  ) +
  theme_plant_scientific() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Violin + Box Plot + Jitter (full combination) ----
ggplot(data_long, aes(x = Group, y = Value)) +
  geom_violin(aes(fill = Group), alpha = 0.3, draw_quantiles = 0.5, size = 0.3) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, alpha = 0.8, size = 0.3) +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.4, color = "grey30") +
  scale_fill_manual(values = okabe_ito) +
  scale_x_discrete(labels = sample_n$label) +
  labs(
    x = "",
    y = "Value (units)",
    fill = "Group"
  ) +
  theme_plant_scientific() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- For small n (<15 per group): beeswarm + box ----
library(ggbeeswarm)
ggplot(data_long, aes(x = Group, y = Value, color = Group)) +
  geom_quasirandom(size = 2.5, alpha = 0.8, width = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, size = 0.3, color = "black") +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    size = 0.3,
    color = "black"
  ) +
  scale_color_manual(values = okabe_ito, guide = "none") +
  labs(
    x = "",
    y = "Value (units)",
    caption = "Crossbar: mean. Error bar: SEM."
  ) +
  theme_plant_scientific()

# ---- Multiple traits in facets ----
ggplot(data_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  facet_wrap(~ Trait, scales = "free_y", ncol = 3) +
  labs(x = "", y = "Value") +
  theme_plant_scientific() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

## Choosing Between Box, Violin, and Beeswarm

| n per group | Recommended | Notes |
|-------------|-------------|-------|
| n < 8 | Beeswarm + mean bar | Show every data point |
| n = 8-30 | Box + jitter | Box shows quartiles; jitter shows individual variation |
| n = 30-100 | Violin + box | Violin shows full distribution shape |
| n > 100 | Violin only or density | Too many points for jitter |

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `draw_quantiles` | c(0.5) | Median line; add c(0.25, 0.75) for quartile lines |
| `width` (boxplot) | 0.6 | Narrower for more groups; wider for 2-3 groups |
| `width` (jitter) | 0.15 | Wider for fewer points; narrower for many groups |
| Jitter alpha | 0.3-0.5 | Lower for large n; higher for small n |

## Plant-Specific Notes

- Phenotype data (plant height, yield, biomass): often right-skewed — consider log-transforming or using violin plot to show the distribution shape
- Multi-environment trials: box plots grouped by environment on x-axis and colored by genotype
- Disease resistance scores: often ordinal (1-9 scale); beeswarm is better than violin for ordinal data
- For time-series (plant growth curves): use line chart, not box plot; see bar-scatter-line.md
- For GWAS phenotype distribution: always include a histogram or box plot of the raw phenotype to show the distribution shape to readers
