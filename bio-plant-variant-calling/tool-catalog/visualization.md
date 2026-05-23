# Variant Calling Visualization

**Goal:** Visual QC of VCF data — variant density, Ti/Tv, missing rate, heterozygosity
**Best for:** All VCF pipelines; essential QC before downstream analysis

## Prerequisites
- R 4.0+
- Packages: data.table, ggplot2, scales, RColorBrewer, gridExtra
- bcftools for stats extraction

## Extract QC Statistics with bcftools

```bash
# Generate comprehensive stats
bcftools stats filtered.vcf.gz > vcf_stats.txt

# Extract per-sample metrics
bcftools stats -s - filtered.vcf.gz > per_sample_stats.txt

# Per-chromosome variant density
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12; do
  bcftools view -r $chr filtered.vcf.gz | bcftools stats | grep "^SN.*number of records" >> chr_counts.txt
done
```

## Variant Density Across Chromosomes

```r
library(data.table)
library(ggplot2)

# Assuming bcftools stats per-chromosome output
# Format: chr, n_variants, chr_length, density_per_kb
density <- fread("chr_density.csv")

ggplot(density, aes(x = reorder(chromosome, -density_per_kb),
                     y = density_per_kb, fill = chromosome)) +
  geom_bar(stat = "identity") +
  labs(title = "Variant Density Across Chromosomes",
       x = "Chromosome", y = "Variants per Kb") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip()

ggsave("outputs/figures/variant_density_chromosomes.pdf", width = 8, height = 6)
```

## Ti/Tv Ratio

```r
# Parse bcftools stats TSTV section
ts_tv_data <- fread("vcf_stats.txt", skip = "TSTV", fill = TRUE)

# Expected Ti/Tv for plants: 2.0-2.5
ggplot(ts_tv_data, aes(x = category, y = ratio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Transition/Transversion Ratio",
       subtitle = "Plant expected range: 2.0-2.5 (dashed lines)",
       x = "Category", y = "Ti/Tv Ratio") +
  theme_minimal()

ggsave("outputs/figures/titv_ratio.pdf", width = 8, height = 5)
```

## Missing Rate Heatmap

```r
library(reshape2)

# Calculate missing rate per variant per sample
# Missing rate = proportion of samples without a call at each site
missing_matrix <- as.matrix(fread("missing_data.csv"))
colnames(missing_matrix) <- paste0("Sample", 1:ncol(missing_matrix))

missing_melt <- melt(missing_matrix)

ggplot(missing_melt, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred",
                      name = "Missing\nRate") +
  labs(title = "Missing Data Rate Heatmap",
       x = "Sample", y = "Chromosome/Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6))

ggsave("outputs/figures/missing_rate_heatmap.pdf", width = 10, height = 8)
```

## Heterozygosity Distribution

```r
# Per-sample heterozygosity
het_data <- fread("per_sample_het.csv")  # columns: sample, het_rate

# Add species type annotation
het_data[, breeding_system := ifelse(species %in% c("rice", "soybean", "wheat"),
                                      "inbred", "outcross")]

ggplot(het_data, aes(x = het_rate, fill = breeding_system)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "gray50") +
  annotate("text", x = 0.05, y = Inf, vjust = 2,
           label = "Inbred threshold (5%)", color = "gray50", size = 3) +
  labs(title = "Heterozygosity Distribution",
       subtitle = "Inbred species: <5% expected; Outcross: 10-30%",
       x = "Heterozygosity Rate", y = "Density") +
  theme_minimal()

ggsave("outputs/figures/heterozygosity_distribution.pdf", width = 8, height = 5)
```

## Multi-Panel Summary Plot

```r
library(gridExtra)

# Combine all QC plots into a single summary figure
p1 <- ggplot(density, aes(x = reorder(chromosome, -density_per_kb),
                           y = density_per_kb)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Variant Density per Kb", x = "", y = "Variants/Kb") +
  theme_minimal() + coord_flip()

p2 <- ggplot(ts_tv_data, aes(x = category, y = ratio)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_hline(yintercept = c(2.0, 2.5), linetype = "dashed", color = "red") +
  labs(title = "Ti/Tv Ratio", x = "", y = "Ratio") +
  theme_minimal()

p3 <- ggplot(missing_melt, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(title = "Missing Rate", x = "Sample", y = "") +
  theme_minimal() + theme(axis.text.y = element_blank())

p4 <- ggplot(het_data, aes(x = het_rate, fill = breeding_system)) +
  geom_density(alpha = 0.6) +
  labs(title = "Heterozygosity", x = "Het Rate", y = "") +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave("outputs/figures/qc_summary.pdf", width = 12, height = 10,
       plot = arrangeGrob(p1, p2, p3, p4, ncol = 2))
```

## Plant-Specific Visualization Notes

- **Polyploid species**: Plot variant density by subgenome (A/B/D for wheat). Subgenome-specific density differences can reveal domestication history.
- **Organelle genomes**: If chloroplast/mitochondrial variants are included, they may dominate density plots due to their small size and high density. Consider separate panels for nuclear vs organelle or normalize by chromosome length.
- **Inbred species**: Near-zero heterozygosity is expected and desirable. If an inbred sample shows >5% heterozygosity, flag it as possible contamination or outcrossing.
- **Pericentromeric regions**: In many plant genomes, variant density drops sharply near centromeres due to low recombination and mapping difficulty. This is expected and should not be flagged as a QC issue.
- **Chromosome naming**: Plant reference chromosomes may be numbered (chr1, chr2, ...) or named (1A, 1B, 1D for wheat). Handle mixed naming gracefully in plots.
