# Ks Histogram + LD Decay Curve

**For**: Comparative genomics (whole-genome duplication events) and population genomics (linkage disequilibrium decay)

## Part 1: Ks Histogram with Peak Fitting

### Data Format Required

| Column | Type | Description |
|--------|------|-------------|
| pair_id | character | Gene pair identifier |
| Ks | numeric | Synonymous substitution rate |
| Ka | numeric | Non-synonymous substitution rate (optional) |
| Ka_Ks | numeric | Ka/Ks ratio (optional) |

### ggplot2 Code

```r
library(ggplot2)
library(dplyr)
library(mclust)  # for Gaussian mixture model peak detection
source("theme/theme_plant_scientific.R")

# Load Ks data
# ks_data <- read.table("ks_values.txt", header = TRUE)

# Filter outliers (Ks > 2 usually indicates saturation)
# ks_data <- ks_data %>% filter(Ks < 2)

# Fit Gaussian mixture model to detect peaks
# ks_values <- ks_data$Ks
# gmm <- Mclust(ks_values, G = 1:5)
# peaks <- gmm$parameters$mean

# ---- Histogram with density overlay ----
ggplot(ks_data, aes(x = Ks)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 200,
    fill = okabe_ito[2],
    alpha = 0.6,
    boundary = 0
  ) +
  geom_density(color = okabe_ito[6], size = 0.8) +
  # Add vertical lines for detected peaks
  # geom_vline(xintercept = peaks, linetype = "dashed", color = "red", size = 0.3) +
  labs(
    x = expression(italic(K)[s] ~ "(synonymous substitutions per site)"),
    y = "Density"
  ) +
  theme_plant_scientific()

# ---- With annotation for WGD events ----
# Label known WGD events if the species has them
wgd_events <- data.frame(
  Ks = c(0.3, 0.8, 1.5),          # example values — replace with real data
  label = c("gamma", "beta", "alpha"),
  color = c(okabe_ito[6], okabe_ito[1], okabe_ito[3])
)

ggplot(ks_data, aes(x = Ks)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 200,
    fill = "grey80",
    alpha = 0.6,
    boundary = 0
  ) +
  geom_density(color = "grey30", size = 0.8) +
  geom_vline(
    data = wgd_events,
    aes(xintercept = Ks, color = label),
    linetype = "dashed",
    size = 0.5
  ) +
  geom_text(
    data = wgd_events,
    aes(x = Ks, y = 0.95 * max(..density..), label = label, color = label),
    angle = 90,
    vjust = -0.5,
    size = 3
  ) +
  scale_color_manual(values = setNames(wgd_events$color, wgd_events$label)) +
  labs(
    x = expression(italic(K)[s]),
    y = "Density"
  ) +
  theme_plant_scientific()
```

## Part 2: LD Decay Curve

### Data Format Required

Two-column format from PLINK `--r2` or PopLDdecay:

| Column | Type | Description |
|--------|------|-------------|
| dist_bp | integer | Physical distance between marker pairs (bp) |
| r2 | numeric | R-squared (LD measure) |
| mean_r2 | numeric | Mean r2 per distance bin (pre-computed or computed in R) |

### ggplot2 Code

```r
library(ggplot2)
library(dplyr)
source("theme/theme_plant_scientific.R")

# Load LD data
# ld_data <- read.table("ld_decay.txt", header = TRUE)

# Bin distances and compute mean r2
ld_binned <- ld_data %>%
  mutate(dist_bin = cut(dist_bp, breaks = seq(0, max(dist_bp) + 1000, by = 1000))) %>%
  group_by(dist_bin) %>%
  summarise(
    mean_dist = mean(dist_bp) / 1000,  # convert to kb
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE)
  )

# ---- LD Decay Curve ----
ggplot(ld_binned, aes(x = mean_dist, y = mean_r2)) +
  geom_line(color = okabe_ito[3], size = 0.8) +
  geom_ribbon(aes(ymin = mean_r2 - sd_r2, ymax = mean_r2 + sd_r2),
              fill = okabe_ito[3], alpha = 0.2) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red", size = 0.4) +
  geom_hline(yintercept = 0.1, linetype = "dotted", color = "grey50", size = 0.3) +
  labs(
    x = "Distance (kb)",
    y = expression(italic(r)^2)
  ) +
  theme_plant_scientific()

# ---- Multiple populations on one plot ----
# ld_data: add a 'Population' column
ggplot(ld_binned, aes(x = mean_dist, y = mean_r2, color = Population)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey50", size = 0.3) +
  scale_color_manual(values = okabe_ito) +
  labs(
    x = "Distance (kb)",
    y = expression(italic(r)^2),
    color = "Population"
  ) +
  theme_plant_scientific()
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| Ks bin width | 0.01 | Narrower for high-quality genomes; wider for fragmented assemblies |
| Ks x-axis max | 2 | Saturation usually occurs >2; adjust for species with recent WGD only |
| LD distance bin | 1 kb | Wider bins for low-density markers; narrower for high-density |
| LD threshold (r2=0.2) | 0.2 | Standard for "useful LD"; adjust to 0.1 for stricter |
| GMM components (K) | 1:5 | For Ks peak detection; more components for complex polyploid history |

## Plant-Specific Notes

### Ks
- Plants have multiple rounds of WGD (paleopolyploidy). Label known events: gamma (eudicot triplication), beta, alpha
- Crop-specific WGD events: maize (~12 MYA), soybean (~13 and ~59 MYA), Brassica (~15-40 MYA), cotton (~60 and ~2 MYA)
- For recently formed polyploids: Ks peaks may overlap; use synteny blocks to separate subgenome contributions
- Ks saturation: monocot nuclear genes saturate ~1.0-1.5, eudicot ~0.8-1.2

### LD Decay
- LD decay distance varies dramatically: maize ~1-2 kb (outcrossing), soybean ~50-150 kb (selfing), rice ~100-200 kb
- Outcrossing species have much faster LD decay than selfing species
- Separate LD curves by chromosome; some chromosomes may have distinct decay patterns
- For GWAS context: mark the average marker density as a vertical line to show whether gaps are covered
