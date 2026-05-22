# GWAS Visualization

**Goal:** Publication-quality Manhattan plot, QQ plot, and significant SNP annotation
**Best for:** All GWAS results, regardless of method used

## Prerequisites
- R 4.0+
- Packages: data.table, ggplot2, qqman, ggrepel

## Manhattan + QQ Plot

```r
library(data.table)
library(ggplot2)
library(qqman)

# Read results (adapt column names based on method)
# GAPIT output: SNP, Chromosome, Position, P.value
# PLINK2 output: #CHROM, POS, P
# BLINK output: similar to GAPIT

gwas <- fread("gwas_results.csv")

# Standardize column names
setnames(gwas, old = c("#CHROM", "CHROM", "Chromosome", "CHR"),
         new = c("CHR", "CHR", "CHR", "CHR"), skip_absent = TRUE)
setnames(gwas, old = c("POS", "BP", "Position"),
         new = c("BP", "BP", "BP"), skip_absent = TRUE)
setnames(gwas, old = c("P", "P.value", "p.value"),
         new = c("P", "P", "P"), skip_absent = TRUE)

gwas <- gwas[!is.na(P), ]

# Calculate genomic inflation factor (λ)
lambda <- median(qchisq(gwas$P, df = 1, lower.tail = FALSE), na.rm = TRUE) /
          qchisq(0.5, df = 1)
cat(sprintf("Genomic inflation λ = %.4f\n", lambda))

# Determine significance threshold
n_snps <- nrow(gwas)
bonf_standard <- 5e-8
bonf_relaxed <- 0.05 / n_snps
threshold <- if (n_snps < 10000) bonf_relaxed else bonf_standard
cat(sprintf("Significance threshold: %.2e\n", threshold))

# Manhattan plot
pdf("outputs/figures/manhattan.pdf", width = 12, height = 6)
manhattan(gwas, chr = "CHR", bp = "BP", p = "P",
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(threshold),
          main = sprintf("Manhattan Plot (λ=%.3f, n=%d SNPs)", lambda, n_snps))
dev.off()

# QQ plot
pdf("outputs/figures/qq_plot.pdf", width = 8, height = 8)
qq(gwas$P, main = sprintf("QQ Plot (λ=%.3f)", lambda))
dev.off()

# Export significant SNPs
sig <- gwas[P < threshold, ][order(P)]
fwrite(sig, "outputs/tables/significant_snps.csv")
cat(sprintf("%d significant SNPs at p < %.2e\n", nrow(sig), threshold))

# Export all results
fwrite(gwas, "outputs/tables/all_snp_results.csv")
```

## Plant-Specific Visualization Notes

- Use species chromosome count for proper x-axis in Manhattan plot (e.g., maize=10, rice=12, wheat=21)
- For polyploid species, color points by subgenome (A/B/D for wheat)
- If the species is unknown, use sequential chromosome numbering
