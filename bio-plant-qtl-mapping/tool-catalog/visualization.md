# QTL Mapping Visualization

**Goal:** Publication-quality LOD profile plots, genetic map plots, QTL interval annotation, and QTL effect plots
**Best for:** All QTL mapping results, regardless of method used (IM, CIM, MQM, ICIM)

## Prerequisites
- R 4.0+
- Packages: qtl, ggplot2, data.table, ggrepel, gridExtra

## LOD Profile Plot

```r
library(qtl)
library(ggplot2)

# --- Method-specific LOD profile extraction ---

# For scanone output (IM / CIM):
lod_data <- as.data.frame(out.cim)
# Columns: chr, pos, lod1, lod2 (for 2-QTL scan)

# Plot LOD profile across chromosomes
plot(out.cim, main = "CIM LOD Profile",
     ylab = "LOD Score", xlab = "Chromosome",
     bandcol = "gray90",
     alternate.chrid = TRUE)

# Add significance threshold
operm <- cim(mycross, ..., n.perm = 1000)  # Permutation results
plot(out.cim)
abline(h = summary(operm)[1], col = "red", lty = 2)  # 5% threshold
abline(h = summary(operm)[2], col = "orange", lty = 3)  # 10% threshold
legend("topright", legend = c("5%", "10%"),
       col = c("red", "orange"), lty = 2:3)

# For stepwiseqtl output: LOD profiles conditional on other QTL
plotLodProfile(refined, main = "MQM QTL LOD Profiles",
               ylab = "Conditional LOD", showallchr = TRUE)
```

## ggplot2 Custom LOD Profile

```r
library(ggplot2)
library(data.table)

lod_dt <- as.data.table(lod_data)

# Add cumulative position for x-axis
chr_info <- lod_dt[, .(max_pos = max(pos)), by = chr]
chr_info[, offset := c(0, cumsum(max_pos[-.N]))]
lod_dt <- merge(lod_dt, chr_info[, .(chr, offset)], by = "chr")
lod_dt[, cum_pos := pos + offset]

# Chromosome midpoints for axis labels
chr_labels <- lod_dt[, .(mid = mean(cum_pos)), by = chr]

# Get threshold from permutation
threshold_5pct <- summary(operm)[1]

# Plot
ggplot(lod_dt, aes(x = cum_pos, y = lod)) +
  geom_rect(aes(xmin = cum_pos, xmax = cum_pos + 1, ymin = 0, ymax = lod,
                fill = chr), alpha = 0.7) +
  geom_hline(yintercept = threshold_5pct, linetype = "dashed",
             color = "red", size = 0.8) +
  scale_x_continuous(breaks = chr_labels$mid, labels = chr_labels$chr) +
  labs(title = "QTL LOD Profile (CIM)",
       x = "Chromosome", y = "LOD Score") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave("outputs/figures/lod_profile.pdf", width = 12, height = 5)
```

## Genetic Map Plot

```r
# Basic genetic map plot
pdf("outputs/figures/genetic_map.pdf", width = 12, height = 8)
plotMap(mycross, main = "Genetic Linkage Map",
        show.marker.names = FALSE, alternate.chrid = TRUE)
dev.off()

# Detailed per-chromosome map
pdf("outputs/figures/genetic_map_detailed.pdf", width = 16, height = 12)
par(mfrow = c(3, 4))  # Adjust based on chromosome count
for (chr_name in names(mycross$geno)) {
  plotMap(mycross, chr = chr_name, show.marker.names = TRUE,
          main = paste("Chromosome", chr_name))
}
dev.off()

# Map with QTL positions annotated
# After QTL analysis
qtl_positions <- data.frame(
  chr = c("1", "3", "5"),
  pos = c(45.5, 78.2, 12.8)
)

plotMap(mycross, show.marker.names = FALSE)
for (i in 1:nrow(qtl_positions)) {
  points(qtl_positions$pos[i], qtl_positions$chr[i],
         pch = 23, bg = "red", cex = 1.5)
}
```

## QTL Effect Plot

```r
# Plot QTL effects: phenotype means by genotype at QTL marker
plotPXG(
  cross = mycross,
  marker = find.marker(mycross, chr = "1", pos = 45.5),
  pheno.col = "trait1",
  main = "QTL Effect: qPH-1-1",
  ylab = "Plant Height (cm)"
)

# For multi-QTL model: effect plot from fitqtl
plot.fitqtl <- function(fit, cross, pheno.col) {
  eff <- fit$ests$ests
  se  <- fit$ests$SEs
  
  # For each QTL, get the nearest marker
  # Extract additive effect and plot
  
  qtl_names <- rownames(eff)
  add_eff <- eff[, "a"]  # Additive effect
  
  # Bar plot of additive effects
  barplot(add_eff, names.arg = qtl_names,
          xlab = "QTL", ylab = "Additive Effect",
          main = paste("QTL Additive Effects for", pheno.col),
          col = ifelse(add_eff > 0, "steelblue", "tomato"))
}

plot.fitqtl(fit_refined, mycross, "trait1")
```

## QTL Interval Annotation on Chromosomes

```r
# Draw chromosomes with QTL confidence intervals marked

draw_qtl_intervals <- function(cross, qtl_summary, filename = "qtl_intervals.pdf") {
  # qtl_summary: data.frame with chr, pos, ci_lo, ci_hi, name columns
  
  pdf(filename, width = 12, height = 8)
  
  # Get chromosome lengths
  chr_lengths <- chrlen(cross)
  n_chr <- length(chr_lengths)
  
  # Plot frame
  plot(0, 0, type = "n", xlim = c(0, max(chr_lengths) + 20),
       ylim = c(0.5, n_chr + 0.5),
       xlab = "Map Position (cM)", ylab = "Chromosome",
       main = "QTL Positions on Genetic Map",
       yaxt = "n")
  
  # Draw chromosomes as lines
  for (i in seq_along(chr_lengths)) {
    lines(c(0, chr_lengths[i]), c(i, i), lwd = 3, col = "gray60")
    axis(2, at = i, labels = names(chr_lengths)[i], las = 2)
  }
  
  # Draw QTL intervals
  for (j in 1:nrow(qtl_summary)) {
    chr_idx <- which(names(chr_lengths) == qtl_summary$chr[j])
    # Confidence interval
    rect(qtl_summary$ci_lo[j], chr_idx - 0.2,
         qtl_summary$ci_hi[j], chr_idx + 0.2,
         col = rgb(1, 0, 0, 0.3), border = "red")
    # Peak position
    points(qtl_summary$pos[j], chr_idx, pch = 19, col = "red", cex = 1.2)
    # Label
    text(qtl_summary$pos[j], chr_idx + 0.35, qtl_summary$name[j],
         cex = 0.7, col = "darkred")
  }
  
  dev.off()
}
```

## Summary Table Export

```r
# Generate QTL summary table
qtl_table <- data.frame(
  QTL_Name = character(),
  Chromosome = character(),
  Peak_cM = numeric(),
  CI_low_cM = numeric(),
  CI_high_cM = numeric(),
  LOD = numeric(),
  PVE_percent = numeric(),
  Additive_Effect = numeric(),
  stringsAsFactors = FALSE
)

# Populate from analysis results
# ...

write.csv(qtl_table, "outputs/tables/qtl_summary.csv", row.names = FALSE)

# Also export for publication (LaTeX)
library(xtable)
xtable_qtl <- xtable(qtl_table,
                     caption = "Identified QTL for Plant Height",
                     label = "tab:qtl_ph")
print(xtable_qtl, file = "outputs/tables/qtl_summary.tex",
      include.rownames = FALSE, booktabs = TRUE)
```

## Plant-Specific Visualization Notes

- Use species-specific chromosome naming (e.g., for rice: chr01-chr12; for wheat: 1A, 1B, 1D, ..., 7D).
- Add physical positions alongside genetic positions when a reference genome is available.
- For multi-environment QTL: use faceted plots (one panel per environment) or color-coded overlay.
- Annotate known QTL or genes from literature for comparison (useful for reviewers).
- For polyploid species: color points by subgenome in Manhattan-style LOD plots.
