# Differential Abundance Analysis (ANCOM-BC)

**Goal:** Identify taxa that are differentially abundant between groups while controlling for compositionality
**Best for:** Comparing microbial communities across plant genotypes, treatments, or compartments

## Prerequisites
- R 4.0+
- R packages: ANCOMBC, phyloseq, qiime2R
- QIIME 2 outputs: feature table (QZA), taxonomy (QZA), metadata (TSV)

## ANCOM-BC Workflow

### Import QIIME 2 Data into R

```r
library(phyloseq)
library(ANCOMBC)
library(qiime2R)

# Read QIIME 2 artifacts
physeq <- qza_to_phyloseq(
  features = "table.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# Inspect data
sample_data(physeq)
tax_table(physeq)[1:5, 1:5]
```

### Run ANCOM-BC at Genus Level

```r
# Aggregate to genus level
physeq_genus <- tax_glom(physeq, taxrank = "Genus")

# Run ANCOM-BC with treatment as grouping variable
out <- ancombc2(
  data = physeq_genus,
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "Treatment + Tissue",
  rand_formula = "(1 | Block)",
  p_adj_method = "holm",
  prv_cut = 0.10,       # Prevalence filter: keep taxa in >10% of samples
  lib_cut = 1000,        # Library size filter
  group = "Treatment",
  struc_zero = TRUE,     # Account for structural zeros
  neg_lb = TRUE,
  alpha = 0.05,
  global = TRUE          # Test overall effect
)

# View significant results
res <- out$res
sig_taxa <- res[res$diff_Treatment == TRUE, ]
head(sig_taxa[, c("taxon", "lfc_Treatment", "p_Treatment", "q_Treatment")])

# Export results
write.csv(sig_taxa, "ancombc_significant_taxa.csv")
```

### Bias-Corrected Abundance Estimation

```r
# Extract bias-corrected abundances
bias_corrected <- out$bias_corrected_log_table
head(bias_corrected)
```

## ALDEx2 (Alternative Method)

```r
library(ALDEx2)

# Prepare count matrix
counts <- as(otu_table(physeq_genus), "matrix")
conditions <- sample_data(physeq_genus)$Treatment

# Run ALDEx2
aldex_out <- aldex(
  reads = t(counts),
  conditions = conditions,
  mc.samples = 128,
  test = "t",
  effect = TRUE,
  denom = "iqlr"  # IQLR denominator for compositionality
)

# Filter significant taxa (effect size > 1 AND p < 0.05)
sig <- aldex_out[abs(aldex_out$effect) > 1 & aldex_out$wi.eBH < 0.05, ]
```

## Method Comparison

| Method | Handles Compositionality | Handles Sparsity | Statistical Framework | Best For |
|--------|-------------------------|------------------|----------------------|----------|
| ANCOM-BC | Yes (bias-corrected) | Yes | Linear model with log-link | Large studies with many groups |
| ALDEx2 | Yes (CLR transform) | Moderate | Bayesian + effect size | Small to medium studies |
| MaAsLin2 | Yes (normalization) | Yes | Mixed effects linear model | Longitudinal studies |
| DESeq2 | No (needs size factors) | No | Negative binomial | Rich communities only |
| LEfSe | No | No | Kruskal-Wallis + LDA | Exploratory only (not recommended as primary) |

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| fix_formula | Fixed effects (treatment groups) |
| rand_formula | Random effects (blocking factors) |
| prv_cut | Prevalence cutoff |
| lib_cut | Minimum library size per sample |
| p_adj_method | Multiple testing correction (holm, BH, BY) |
| struc_zero | Model structural zeros separately |
| group | Primary group variable for comparison |
| alpha | Significance level |

## Plant-Specific Considerations

- Plant genotype as random effect: include genotype/bloc in `rand_formula` when samples are from multiple cultivars
- Compartment effect: root compartment (rhizosphere vs endosphere) is typically the dominant source of variation; include it as a fixed effect
- Soil batch: if soils were collected at different sites/times, include site as a random effect
- Time series: for plant development studies, use MaAsLin2 with time as a continuous variable instead of ANCOM-BC

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No significant taxa" | Low effect size or small sample size | Relax prv_cut; check statistical power |
| "Design matrix singular" | Confounded metadata columns | Simplify formula; check for collinearity |
| "Error in ancombc2: negative values" | Non-count data in assay | Ensure feature table contains raw counts, not relative abundance |
| All taxa significant | Multiple testing issue | Use stricter p_adj_method (holm vs BH) |
