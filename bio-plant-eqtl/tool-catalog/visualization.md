# eQTL Visualization

**Goal:** Generate publication-quality figures for eQTL results — Manhattan plots,
cis/trans classification summaries, colocalization plots, and multi-tissue sharing
visualizations.

**Best for:** Summarizing eQTL analysis results for reports, presentations, and
publications.

## Prerequisites

- R 4.0+, ggplot2, qqman (for Manhattan), ComplexHeatmap or pheatmap
- eQTL results from MatrixEQTL or similar tool
- Gene and SNP annotations

## Workflow

### 1. eQTL Manhattan Plot (per Chromosome or Genome-Wide)

```r
library(ggplot2)
library(dplyr)

# Prepare data: SNP positions, -log10(p), chromosome
eqtl_plot_data <- cis_eqtl %>%
  mutate(
    log10p = -log10(pvalue),
    sig = ifelse(fdr < 0.05, "FDR < 0.05",
          ifelse(pvalue < 1e-5, "Suggestive", "NS"))
  )

# Per-chromosome Manhattan (for a single eGene)
plot_eqtl_manhattan <- function(eqtl_data, gene_name, cis_window = 1e6) {
  gene_info <- gene_loc[gene_loc$geneid == gene_name, ]
  x_limits <- c(gene_info$start - cis_window, gene_info$start + cis_window)

  ggplot(eqtl_data, aes(x = pos, y = log10p, color = sig)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(
      values = c("FDR < 0.05" = "#B2182B",
                 "Suggestive" = "#2166AC",
                 "NS" = "grey70")
    ) +
    geom_vline(xintercept = gene_info$start, linetype = "dashed",
               color = "darkgreen", alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted",
               color = "black") +
    annotate("text", x = gene_info$start, y = max(eqtl_data$log10p),
             label = gene_name, vjust = -1, color = "darkgreen",
             fontface = "italic") +
    labs(x = paste0("Position on ", gene_info$chr),
         y = expression(-log[10](italic(P))),
         title = paste("cis-eQTL for", gene_name)) +
    theme_bw() +
    theme(legend.position = "bottom")
}

# Example
plot_eqtl_manhattan(
  subset(eqtl_plot_data, gene == "Os01g0100100" & chr == "Chr1"),
  "Os01g0100100"
)
```

### 2. Cis vs Trans Summary

```r
# Bar chart: number of cis vs trans eQTLs
eqtl_summary <- eqtl %>%
  mutate(type = factor(type, levels = c("cis", "trans"))) %>%
  group_by(type) %>%
  summarise(
    n_total = n(),
    n_significant = sum(fdr < 0.05, na.rm = TRUE)
  )

ggplot(eqtl_summary, aes(x = type, y = n_total, fill = type)) +
  geom_col(alpha = 0.8) +
  geom_col(aes(y = n_significant), fill = "#B2182B", alpha = 0.9, width = 0.5) +
  geom_text(aes(label = paste0(n_significant, " sig")),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("cis" = "#2166AC", "trans" = "#4393C3")) +
  labs(x = "eQTL Type", y = "Number of SNP-Gene Pairs",
       title = "Cis vs Trans eQTL Classification",
       subtitle = paste0("Significant at FDR < 0.05")) +
  theme_bw()
```

### 3. eGene Distribution Across Chromosomes

```r
# Count eGenes per chromosome
egene_chr <- egenes %>%
  left_join(gene_loc, by = c("gene" = "geneid")) %>%
  count(chr, name = "n_eGenes")

# Normalize by total genes per chromosome
gene_count_per_chr <- gene_loc %>% count(chr, name = "total_genes")
egene_chr <- egene_chr %>%
  left_join(gene_count_per_chr, by = "chr") %>%
  mutate(proportion = n_eGenes / total_genes)

ggplot(egene_chr, aes(x = reorder(chr, as.numeric(gsub("Chr", "", chr))),
                      y = proportion)) +
  geom_col(fill = "#2166AC", alpha = 0.8) +
  geom_text(aes(label = n_eGenes), vjust = -0.5, size = 3) +
  labs(x = "Chromosome", y = "Proportion of Genes that are eGenes",
       title = "eGene Density per Chromosome") +
  theme_bw()
```

### 4. Colocalization Posterior Probability Plot

```r
# Plot coloc results for a single locus
plot_coloc_results <- function(coloc_summary, locus_name) {
  coloc_df <- data.frame(
    Hypothesis = c("H0: Neither", "H1: eQTL only",
                   "H2: GWAS only", "H3: Independent",
                   "H4: Shared"),
    Posterior = as.numeric(coloc_summary[c("PP.H0.abf", "PP.H1.abf",
                                           "PP.H2.abf", "PP.H3.abf",
                                           "PP.H4.abf"), ])
  )

  ggplot(coloc_df, aes(x = Hypothesis, y = Posterior, fill = Hypothesis)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "#B2182B") +
    annotate("text", x = 5, y = 0.78, label = "PP.H4 > 0.75 threshold",
             color = "#B2182B", hjust = 1, size = 3.5) +
    scale_fill_manual(values = c("grey80", "#4393C3", "#F4A582",
                                 "#92C5DE", "#2166AC")) +
    labs(title = paste("Colocalization:", locus_name),
         y = "Posterior Probability") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
}
```

### 5. Multi-Tissue eQTL Sharing Heatmap

```r
library(pheatmap)

# Correlation of eQTL effect sizes across tissues
# (from mashR output or pairwise comparisons)
tissue_cor <- cor(eqtl_effects_matrix, method = "spearman")

pheatmap(tissue_cor,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("#B2182B", "white", "#2166AC"))(100),
         main = "eQTL Effect Size Correlation Across Tissues",
         fontsize = 11,
         border_color = NA)
```

### 6. QQ Plot for eQTL P-values

```r
plot_eqtl_qq <- function(pvalues, title = "eQTL QQ Plot") {
  n <- length(pvalues)
  observed <- sort(pvalues)
  expected <- (1:n) / (n + 1)

  # Genomic inflation factor (lambda)
  chisq <- qchisq(1 - observed, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)

  df <- data.frame(
    expected_log10 = -log10(expected),
    observed_log10 = -log10(observed)
  )

  ggplot(df, aes(x = expected_log10, y = observed_log10)) +
    geom_point(alpha = 0.3, size = 0.5, color = "#2166AC") +
    geom_abline(slope = 1, intercept = 0, color = "#B2182B", linetype = "dashed") +
    annotate("text", x = 0.5, y = max(df$observed_log10) * 0.95,
             label = paste("lambda =", round(lambda, 3)),
             hjust = 0, size = 4) +
    labs(x = expression(Expected ~ -log[10](italic(P))),
         y = expression(Observed ~ -log[10](italic(P))),
         title = title) +
    theme_bw()
}

# Apply to cis-eQTL p-values
plot_eqtl_qq(cis_eqtl$pvalue, "cis-eQTL QQ Plot")
```

### 7. Effect Size vs MAF

```r
# Identify high-effect, rare-variant eQTLs
ggplot(cis_eqtl, aes(x = maf, y = abs(beta), color = fdr < 0.05)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey70"),
                     labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS")) +
  scale_x_log10() +
  labs(x = "Minor Allele Frequency (log10 scale)",
       y = "|Effect Size (beta)|",
       title = "eQTL Effect Size vs MAF",
       color = "") +
  theme_bw() +
  theme(legend.position = "bottom")
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| Manhattan point size | 0.5-1.0 | Balance visibility and overplotting |
| Manhattan alpha | 0.3-0.6 | Helps see point density |
| Heatmap clustering method | "complete" or "ward.D2" | Standard for gene expression patterns |
| QQ plot lambda annotation | median-based | Standard genomic inflation metric |
| Coloc threshold line | PP.H4 = 0.75 | Standard colocalization threshold |
| FDR significance line | -log10(FDR 0.05 equivalent) | Context-dependent |

## Plant-Specific Notes

- **Chromosome naming**: Plant genomes use diverse naming (Chr01, chr1A, 1, NC_001).
  Normalize chromosome names before plotting for consistent sort order.
- **Large genomes**: Wheat and barley have chromosome-level data with many scaffolds.
  Consider plotting only assembled chromosomes and labeling unplaced scaffolds as "Un".
- **Polyploid subgenome separation**: When plotting eQTLs for a polyploid, use
  faceted Manhattan plots by subgenome, or use distinct color palettes per subgenome.
- **Multi-panel figures**: Plant journals often require combining Manhattan, QQ,
  and table in a single figure. Use `cowplot::plot_grid()` or `patchwork`.
- **Gene naming conventions**: Italicize gene names in plot labels following plant
  nomenclature. Note that different species have different conventions (Os for rice,
  AT for Arabidopsis, Zm for maize, Ta for wheat).
- **Environmental context**: When plotting eQTL from multiple environments, always
  show environment labels clearly. eQTL Manhattan plots without environment context
  are difficult to interpret.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Manhattan plot: chromosomes sorted as 1,10,11,...,2,20 | String sorting of chr names | Use `mixedsort` or strip non-numeric characters |
| Heatmap: all correlations ~0 | Different eQTLs in different tissues (biological) | Check if common eQTLs exist; plot eQTL overlap count |
| QQ plot: extreme tail inflation | Missing covariates or population structure | Add more PEER factors; check for batch effects |
| Coloc plot: H3 dominates | Independent signals (common in plants with long LD) | Try conditional analysis or adjust window |
| Too many points in scatter | Large eQTL dataset | Use `geom_hex()` or `geom_bin2d()` for density |
