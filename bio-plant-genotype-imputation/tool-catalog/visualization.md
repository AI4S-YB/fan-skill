# Genotype Imputation Visualization

**Goal:** Generate publication-quality plots for imputation accuracy and data summary
**Best for:** All imputation results, regardless of tool used

## Prerequisites

- R 4.0+
- Packages: data.table, ggplot2, gridExtra, scales, RColorBrewer

## Imputation Accuracy Distribution

```r
library(data.table)
library(ggplot2)

# Read variant metrics (output from post-imputation QC)
# Columns: CHR, POS, AF, DR2 (or R2, or INFO)

metrics <- fread("variant_metrics.txt",
                  col.names = c("CHR", "POS", "AF", "DR2"))
metrics$MAF <- pmin(metrics$AF, 1 - metrics$AF)
metrics <- metrics[!is.na(DR2), ]

# 1. DR2 Distribution Histogram
pdf("outputs/figures/dr2_distribution.pdf", width = 8, height = 6)
ggplot(metrics, aes(x = DR2)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(0.3, 0.5, 0.8), linetype = "dashed",
             color = c("red", "orange", "darkgreen")) +
  annotate("text", x = c(0.3, 0.5, 0.8), y = Inf, vjust = 2,
           label = c("R²=0.3", "R²=0.5", "R²=0.8"),
           color = c("red", "orange", "darkgreen")) +
  labs(title = "Imputation Accuracy (DR2) Distribution",
       subtitle = sprintf("Mean DR2 = %.4f | %d variants",
                          mean(metrics$DR2, na.rm = TRUE), nrow(metrics)),
       x = "Dosage R² (DR2)", y = "Count") +
  theme_bw()
dev.off()
```

## DR2 by MAF Bin

```r
# 2. DR2 vs MAF scatter (subsample for large datasets)
set.seed(42)
plot_data <- metrics[sample(.N, min(.N, 50000)), ]

pdf("outputs/figures/dr2_by_maf.pdf", width = 8, height = 6)
ggplot(plot_data, aes(x = MAF, y = DR2)) +
  geom_point(alpha = 0.1, size = 0.5, color = "steelblue") +
  geom_smooth(method = "gam", color = "red") +
  labs(title = "Imputation Accuracy vs Minor Allele Frequency",
       subtitle = sprintf("n = %d variants (subsampled if > 50K)", nrow(plot_data)),
       x = "Minor Allele Frequency (MAF)", y = "Dosage R² (DR2)") +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme_bw()
dev.off()

# 3. Boxplot by MAF bin
metrics$MAF_bin <- cut(metrics$MAF,
                        breaks = c(0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5),
                        labels = c("<0.5%", "0.5-1%", "1-5%", "5-10%", "10-20%", ">20%"))

pdf("outputs/figures/dr2_boxplot_by_maf_bin.pdf", width = 10, height = 6)
ggplot(metrics, aes(x = MAF_bin, y = DR2, fill = MAF_bin)) +
  geom_boxplot(outlier.alpha = 0.1) +
  labs(title = "DR2 by MAF Bin",
       x = "MAF Bin", y = "Dosage R² (DR2)") +
  theme_bw() + theme(legend.position = "none")
dev.off()
```

## Before/After Marker Density Comparison

```r
# 4. Marker density comparison: pre vs post imputation
pre_snps <- fread("pre_imputation_snp_positions.txt",
                   col.names = c("CHR", "POS"))
post_snps <- fread("variant_metrics.txt",
                    col.names = c("CHR", "POS", "AF", "DR2"))

# Count SNPs per 1 Mb window
count_snps <- function(data, chr_col = "CHR", pos_col = "POS") {
  data$window <- floor(data[[pos_col]] / 1e6)
  as.data.table(table(data[[chr_col]], data$window))
}

pre_counts <- count_snps(pre_snps)
post_counts <- count_snps(post_snps)

pdf("outputs/figures/snp_density_comparison.pdf", width = 14, height = 8)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
barplot(pre_counts$N, main = "Pre-Imputation SNP Density (per 1Mb)",
        xlab = "Genomic Window (1Mb)", ylab = "SNP Count",
        col = "gray60", border = NA)
barplot(post_counts$N, main = "Post-Imputation SNP Density (per 1Mb)",
        xlab = "Genomic Window (1Mb)", ylab = "SNP Count",
        col = "steelblue", border = NA)
dev.off()
```

## Per-Chromosome DR2 Summary

```r
# 5. DR2 per chromosome
chr_summary <- metrics[, .(
  mean_DR2 = mean(DR2, na.rm = TRUE),
  sd_DR2 = sd(DR2, na.rm = TRUE),
  n_variants = .N,
  prop_dr2_gt_0.3 = sum(DR2 > 0.3, na.rm = TRUE) / .N,
  prop_dr2_gt_0.5 = sum(DR2 > 0.5, na.rm = TRUE) / .N,
  prop_dr2_gt_0.8 = sum(DR2 > 0.8, na.rm = TRUE) / .N
), by = CHR]

pdf("outputs/figures/dr2_per_chromosome.pdf", width = 12, height = 6)
ggplot(chr_summary, aes(x = factor(CHR), y = mean_DR2)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = c(0.3, 0.5, 0.8), linetype = "dashed",
             color = c("red", "orange", "darkgreen")) +
  geom_errorbar(aes(ymin = mean_DR2 - sd_DR2, ymax = mean_DR2 + sd_DR2),
                width = 0.2) +
  labs(title = "Mean Imputation Accuracy (DR2) per Chromosome",
       subtitle = "Error bars: +/- 1 SD",
       x = "Chromosome", y = "Mean DR2") +
  theme_bw()
dev.off()
```

## Pre/Post Summary Table

```r
# 6. Summary statistics table
summary_stats <- data.frame(
  Metric = c("Pre-imputation variants",
             "Post-imputation variants",
             "Fold increase",
             "Mean DR2",
             "Median DR2",
             "% variants DR2 > 0.3",
             "% variants DR2 > 0.5",
             "% variants DR2 > 0.8",
             "Samples"),
  Value = c(
    nrow(pre_snps),
    nrow(metrics),
    round(nrow(metrics) / nrow(pre_snps), 1),
    round(mean(metrics$DR2, na.rm = TRUE), 4),
    round(median(metrics$DR2, na.rm = TRUE), 4),
    round(sum(metrics$DR2 > 0.3, na.rm = TRUE) / nrow(metrics) * 100, 1),
    round(sum(metrics$DR2 > 0.5, na.rm = TRUE) / nrow(metrics) * 100, 1),
    round(sum(metrics$DR2 > 0.8, na.rm = TRUE) / nrow(metrics) * 100, 1),
    "N samples from input"
  )
)

fwrite(summary_stats, "outputs/tables/imputation_summary.csv")

# Print to console
print(summary_stats)
```

## Plant-Specific Visualization Notes

- **Polyploid species**: If subgenome information is available, facet by subgenome in DR2 plots
- **Multi-chromosome species**: Use species-specific chromosome count and naming
- **Large genomes**: Adjust window size for density plots (larger windows for larger genomes)
- **Centromeric regions**: Overlay centromere positions on Manhattan-like DR2 plots to flag unreliable regions
- **Color scheme**: Use `scale_fill_brewer(palette = "Set2")` for multi-category plots (subgenomes, populations)

## Missing Rate Heatmap (Optional)

```r
# 7. Missing rate heatmap (if pre-imputation missing data is available)
library(reshape2)

# Assuming a matrix: rows = samples, cols = genomic bins, values = missing rate
# missing_matrix <- ...

# pdf("outputs/figures/missing_rate_heatmap.pdf", width = 12, height = 8)
# heatmap.2(missing_matrix, trace = "none", col = colorRampPalette(c("white", "red"))(100),
#           main = "Per-Sample Missing Rate by Genomic Region",
#           xlab = "Genomic Region", ylab = "Sample",
#           margins = c(10, 5))
# dev.off()
```
