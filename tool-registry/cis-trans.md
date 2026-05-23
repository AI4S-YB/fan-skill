# Cis-Trans eQTL Classification

**Goal:** Distinguish cis-acting from trans-acting eQTLs based on genomic distance
between SNP and target gene.

**Best for:** Classifying eQTL results after MatrixEQTL or similar tools, and applying
appropriate FDR thresholds for each class.

## Prerequisites

- Completed eQTL mapping (MatrixEQTL output or equivalent)
- SNP location file (SNP ID, chromosome, position)
- Gene location file (gene ID, chromosome, TSS)

## Workflow

### Step 1: Standard Cis/Trans Classification

```r
library(dplyr)

# Load eQTL results and location files
eqtl_results <- read.table("cis_eqtl_results.txt", header = TRUE,
                           stringsAsFactors = FALSE)
snp_loc <- read.table("snp_location.txt", header = TRUE)
gene_loc <- read.table("gene_location.txt", header = TRUE)

# Merge with location info
eqtl <- eqtl_results %>%
  left_join(snp_loc, by = c("snps" = "snpid")) %>%
  rename(snp_chr = chr, snp_pos = pos) %>%
  left_join(gene_loc, by = c("gene" = "geneid")) %>%
  rename(gene_chr = chr, gene_tss = start)

# Classify: cis = same chr AND distance <= 1 Mb
eqtl$type <- with(eqtl, ifelse(
  snp_chr == gene_chr & abs(snp_pos - gene_tss) <= 1e6,
  "cis", "trans"
))

# Summary
cat("cis-eQTLs:", sum(eqtl$type == "cis"), "\n")
cat("trans-eQTLs:", sum(eqtl$type == "trans"), "\n")
```

### Step 2: Per-Gene FDR Correction (eGene Calling)

```r
# For each gene, keep the most significant cis-SNP
top_cis_per_gene <- eqtl %>%
  filter(type == "cis") %>%
  group_by(gene) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

# BH FDR correction across genes
top_cis_per_gene$fdr <- p.adjust(top_cis_per_gene$pvalue, method = "BH")

# eGenes: genes with at least one significant cis-eQTL
egenes <- top_cis_per_gene %>%
  filter(fdr < 0.05)

cat("cis-eGenes (FDR < 0.05):", nrow(egenes), "\n")
```

### Step 3: Trans-eQTL FDR (Genome-Wide)

```r
# Trans correction is more stringent (genome-wide scope)
trans_eqtl <- eqtl %>% filter(type == "trans")

# Bonferroni correction for trans
n_trans_tests <- nrow(trans_eqtl)
trans_eqtl$bonferroni <- trans_eqtl$pvalue * n_trans_tests

# Or FDR for trans
trans_eqtl$fdr <- p.adjust(trans_eqtl$pvalue, method = "BH")

trans_significant <- trans_eqtl %>%
  filter(fdr < 0.05)

cat("Significant trans-eQTLs (FDR < 0.05):", nrow(trans_significant), "\n")
```

### Step 4: Distance-Based Visualization

```r
library(ggplot2)

# Volcano plot: effect size vs significance, colored by type
ggplot(eqtl, aes(x = beta, y = -log10(pvalue), color = type)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_manual(values = c("cis" = "#2166AC", "trans" = "#B2182B")) +
  geom_hline(yintercept = -log10(0.05 / nrow(eqtl)), linetype = "dashed") +
  labs(x = "Effect Size (beta)", y = "-log10(P-value)",
       title = "Cis vs Trans eQTL") +
  theme_bw()
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| cis distance threshold | 1 Mb | Standard for most plant species |
| cis distance (small genome) | 500 kb | Arabidopsis (~135 Mb), rice (~400 Mb) |
| cis distance (large genome) | 1-2 Mb | Wheat (~17 Gb), barley (~5 Gb) |
| cis FDR threshold | 0.05 | Standard eGene definition |
| trans FDR threshold | 0.05 | More stringent; consider 0.01 for hotspot claims |
| min eGene distance to TSS | < 10 kb median | Good data quality indicator |

## Plant-Specific Notes

- **Self-pollinating species** (rice, soybean, wheat): LD extends farther than
  outcrossers. cis window of 1 Mb is conservative but appropriate.
- **Cross-pollinating species** (maize, sunflower): LD decays more rapidly.
  500 kb cis window may be sufficient.
- **Centromeric regions**: LD blocks can span >10 Mb. SNPs in these regions may
  be misclassified as trans. Flag centromeric eQTLs separately.
- **Polyploid homeolog confusion**: A cis-eQTL in subgenome A may appear as a
  trans-eQTL in subgenome B due to cross-mapping. Filter with MAPQ >= 30.
- **Tissue-specific cis/trans ratios**: Root tissues often show higher trans-eQTL
  ratios than leaf tissues — this is biological, not a QC issue.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Very low trans-eQTL count | Low sample size (n < 100) | Expected; trans requires n >= 200 |
| Extremely high trans count | Cross-homeolog mapping in polyploids | Filter SNPs by subgenome mapping quality |
| Cis-eQTL distance distribution bimodal | Centromeric LD artifacts | Remove centromeric SNPs or flag them |
| FDR-adjusted p-values all > 0.05 | Low power or over-correction | Check n; consider suggestive threshold (FDR < 0.1) |
| Negative beta for all cis-eQTL | Reference allele mis-assigned | Check allele coding direction in genotype data |
