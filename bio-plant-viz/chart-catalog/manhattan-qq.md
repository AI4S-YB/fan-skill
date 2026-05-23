# Manhattan Plot + QQ Plot

**For**: GWAS results visualization

## Data Format Required

Two input files or data frames:

### Manhattan plot data
| Column | Type | Description |
|--------|------|-------------|
| SNP | character | SNP/marker ID |
| CHR | integer | Chromosome number |
| BP | integer | Base-pair position |
| P | numeric | p-value from association test |

### QQ plot data
Uses the same p-value column. The qqman package computes expected p-values internally.

## ggplot2 Code

```r
library(qqman)
library(ggplot2)
library(patchwork)
source("theme/theme_plant_scientific.R")

# Load GWAS results
# gwas <- read.table("gwas_results.txt", header = TRUE)

# ---- Manhattan Plot ----
# Using qqman::manhattan() for speed, or ggplot2 for full control

# Option 1: qqman (fast, less customizable)
png("manhattan_qqman.png", width = 10, height = 4, units = "in", res = 300)
manhattan(gwas, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          col = okabe_ito[1:2],  # alternating chromosome colors
          main = "")
dev.off()

# Option 2: ggplot2 (fully customizable, recommended for publication)
# Prepare data: compute cumulative position
gwas <- gwas[order(gwas$CHR, gwas$BP), ]
chr_offsets <- gwas %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(chr_offset = lag(cumsum(as.numeric(chr_len)), default = 0))

gwas <- gwas %>%
  left_join(chr_offsets, by = "CHR") %>%
  mutate(pos_cum = BP + chr_offset)

# Chromosome label positions
axis_pos <- chr_offsets %>%
  mutate(midpoint = chr_offset + chr_len / 2)

# Significance thresholds
suggestive <- -log10(1e-5)
genomewide  <- -log10(5e-8)

# Color assignment
gwas$chr_color <- factor(gwas$CHR %% 2)

# Top SNPs to highlight (optional)
top_snps <- gwas %>% filter(P < 5e-8) %>% slice_min(P, n = 10)

ggplot(gwas, aes(x = pos_cum, y = -log10(P))) +
  geom_point(aes(color = chr_color), size = 0.5, alpha = 0.7) +
  scale_color_manual(values = c(okabe_ito[1], okabe_ito[2]), guide = "none") +
  geom_hline(yintercept = suggestive, linetype = "dashed", color = "grey40", size = 0.3) +
  geom_hline(yintercept = genomewide, linetype = "dashed", color = "red", size = 0.4) +
  scale_x_continuous(breaks = axis_pos$midpoint, labels = axis_pos$CHR) +
  geom_text_repel(data = top_snps, aes(label = SNP), size = 2.5, max.overlaps = 15) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme_plant_scientific() +
  theme(aspect.ratio = 0.3)

# ---- QQ Plot ----
# Option 1: qqman
qq(gwas$P, main = "")

# Option 2: ggplot2
gwas <- gwas %>%
  arrange(P) %>%
  mutate(
    observed = -log10(P),
    expected = -log10(ppoints(n()))
  )

ggplot(gwas, aes(x = expected, y = observed)) +
  geom_abline(slope = 1, intercept = 0, color = "red", size = 0.5) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(x = expression(Expected ~ -log[10](italic(P))),
       y = expression(Observed ~ -log[10](italic(P)))) +
  theme_plant_scientific() +
  coord_fixed()

# ---- Combined Figure (patchwork) ----
# manhattan_plot + qq_plot + plot_layout(heights = c(2, 1), ncol = 1)
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `suggestiveline` | 1e-5 | Suggestive significance threshold; adjust for marker density |
| `genomewideline` | 5e-8 | Bonferroni-corrected genome-wide significance |
| Point size | 0.5 | Smaller for many SNPs; larger for few markers |
| Point alpha | 0.5 | Lower for dense; higher for sparse |

## Plant-Specific Notes

- For polyploid species (wheat, cotton, canola): plot subgenomes with different colors or separate panels labeled A/B/D
- For species with many small chromosomes: consider reducing point size and alpha significantly
- If GWAS was run with multiple models (MLM, FarmCPU, BLINK): use different panels, not overlay on same plot
- LD decay distance: add vertical bands showing LD decay range around significant peaks
- For GAPIT/FarmCPU output: column names may differ; check and rename as needed
