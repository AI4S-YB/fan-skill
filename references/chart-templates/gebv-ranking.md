# GEBV Ranking Bar Chart

**For**: Genomic selection results — ranking individuals by genomic estimated breeding value

## Data Format Required

| Column | Type | Description |
|--------|------|-------------|
| ID | character | Genotype/line identifier |
| GEBV | numeric | Genomic estimated breeding value |
| SE_GEBV | numeric | Standard error of GEBV (optional) |
| Reliability | numeric | Prediction reliability 0-1 (optional) |
| Population | character/factor | Training/validation set or subpopulation |
| Phenotype | numeric | Observed phenotype (optional, for comparison) |

## ggplot2 Code

```r
library(ggplot2)
library(dplyr)
source("theme/theme_plant_scientific.R")

# Load GEBV data
# gebv <- read.table("gebv_results.txt", header = TRUE)

# Sort by GEBV descending
gebv <- gebv %>%
  arrange(desc(GEBV)) %>%
  mutate(
    rank = 1:n(),
    top = rank <= 10  # highlight top 10
  )

# ---- GEBV Ranking Bar Chart ----
ggplot(gebv, aes(x = reorder(ID, GEBV), y = GEBV)) +
  geom_col(aes(fill = top), width = 0.7) +
  geom_errorbar(
    aes(ymin = GEBV - SE_GEBV, ymax = GEBV + SE_GEBV),
    width = 0.3,
    size = 0.3
  ) +
  scale_fill_manual(
    values = c("TRUE" = okabe_ito[6], "FALSE" = "grey70"),
    guide = "none"
  ) +
  labs(
    x = "Genotype",
    y = "GEBV",
    caption = "Error bars: SE of GEBV"
  ) +
  theme_plant_scientific() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
  )

# ---- GEBV Distribution (density/histogram) ----
ggplot(gebv, aes(x = GEBV)) +
  geom_density(fill = okabe_ito[2], alpha = 0.4, color = okabe_ito[2], size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.4) +
  labs(
    x = "GEBV",
    y = "Density"
  ) +
  theme_plant_scientific()

# ---- With reliability coloring ----
ggplot(gebv, aes(x = rank, y = GEBV)) +
  geom_col(aes(fill = Reliability), width = 0.7) +
  scale_fill_viridis_c(option = "D", direction = 1,
                       name = "Reliability") +
  labs(
    x = "Genotype Rank",
    y = "GEBV"
  ) +
  theme_plant_scientific()

# ---- Prediction vs Observation Scatter (for validation) ----
# gebv: must have both GEBV and Phenotype columns
cor_val <- cor(gebv$GEBV, gebv$Phenotype, use = "complete.obs")
r2_val <- cor_val^2

ggplot(gebv, aes(x = Phenotype, y = GEBV)) +
  geom_point(size = 2, alpha = 0.7, color = okabe_ito[3]) +
  geom_smooth(method = "lm", se = TRUE, color = okabe_ito[6], size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  annotate(
    "text",
    x = min(gebv$Phenotype, na.rm = TRUE) + 0.05 * diff(range(gebv$Phenotype, na.rm = TRUE)),
    y = max(gebv$GEBV, na.rm = TRUE) - 0.05 * diff(range(gebv$GEBV, na.rm = TRUE)),
    label = paste0("r = ", round(cor_val, 3), "\nR^2 = ", round(r2_val, 3)),
    hjust = 0,
    size = 3
  ) +
  labs(
    x = "Observed Phenotype",
    y = "Predicted GEBV"
  ) +
  theme_plant_scientific() +
  coord_fixed()
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `top` threshold | 10 | Number of top genotypes to highlight |
| Error bar type | SE | Always specify in caption; PEV (prediction error variance) is also common |
| `coord_fixed()` | TRUE for scatter | Keeps 1:1 aspect ratio for prediction-observation comparison |
| Correlation metric | Pearson r | For prediction accuracy; use Spearman for rank-based evaluation |

## Plant-Specific Notes

- Selection intensity: add a horizontal line at the truncation point (e.g., select top 20%)
- Breeding cycle: if comparing multiple cycles/years, use grouped bars or facets
- Multi-trait selection: show GEBV ranking panels side by side for different traits
- Multi-environment: use box plot or confidence interval overlap to show GxE
- For hybrid breeding: separate male and female GEBV; plot general combining ability (GCA) instead
- Reliability distribution: include a reliability histogram as an inset or adjacent panel
- If comparing multiple GS models: use a grouped bar chart with model on x-axis and correlation on y-axis
