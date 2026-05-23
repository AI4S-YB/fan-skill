# Colocalization Analysis (coloc)

**Goal:** Determine whether eQTL and GWAS signals at the same locus share a causal
variant using Bayesian colocalization.

**Best for:** Integrating eQTL results with existing GWAS summary statistics to
prioritize candidate causal genes.

## Prerequisites

- R 4.0+, coloc package
- eQTL summary statistics (SNP ID, beta, SE, p-value, MAF, sample size)
- GWAS summary statistics (same format)
- SNP location (chromosome, position)
- LD reference or individual-level genotype data (for sensitivity analysis)

## Basic Usage

```r
library(coloc)

# --- Prepare eQTL dataset ---
# Required columns: beta, varbeta (SE^2), snp, position, MAF, N, type="quant"
eqtl_data <- list(
  beta = eqtl_sumstats$beta,
  varbeta = eqtl_sumstats$se^2,
  snp = eqtl_sumstats$snp_id,
  position = eqtl_sumstats$position,
  type = "quant",           # expression is quantitative
  MAF = eqtl_sumstats$maf,
  N = eqtl_sumstats$n_samples,
  sdY = 1                   # standardized expression
)

# --- Prepare GWAS dataset ---
gwas_data <- list(
  beta = gwas_sumstats$beta,
  varbeta = gwas_sumstats$se^2,
  snp = gwas_sumstats$snp_id,
  position = gwas_sumstats$position,
  type = "quant",           # or "cc" for case-control
  MAF = gwas_sumstats$maf,
  N = gwas_sumstats$n_samples,
  sdY = 1                   # standardized phenotype
)

# --- Run coloc ---
# Check dataset alignment
check_dataset(eqtl_data)
check_dataset(gwas_data)

# Colocalization analysis
coloc_result <- coloc.abf(
  dataset1 = eqtl_data,
  dataset2 = gwas_data,
  p1 = 1e-4,   # prior: SNP associated with expression
  p2 = 1e-4,   # prior: SNP associated with trait
  p12 = 1e-5   # prior: SNP associated with both
)

# Results summary
coloc_result$summary
# H0: neither trait has a genetic association
# H1: only trait 1 (eQTL) has an association
# H2: only trait 2 (GWAS) has an association
# H3: both traits, different causal variants
# H4: both traits, shared causal variant

# Posterior probabilities
print(coloc_result$summary)
```

## Batch Colocalization for Multiple Loci

```r
# For each GWAS locus with overlapping cis-eQTL signal:

run_coloc_locus <- function(locus_chr, locus_start, locus_end,
                            eqtl_full, gwas_full, window = 500e3) {

  # Extract variants in window around locus
  eqtl_region <- subset(eqtl_full,
    chr == locus_chr &
    pos >= locus_start - window &
    pos <= locus_end + window)

  gwas_region <- subset(gwas_full,
    chr == locus_chr &
    pos >= locus_start - window &
    pos <= locus_end + window)

  # Intersect by SNP ID
  shared_snps <- intersect(eqtl_region$snp_id, gwas_region$snp_id)

  if (length(shared_snps) < 20) {
    warning("Locus ", locus_chr, ":", locus_start, " has <20 shared SNPs")
    return(NULL)
  }

  eqtl_region <- eqtl_region[eqtl_region$snp_id %in% shared_snps, ]
  gwas_region <- gwas_region[gwas_region$snp_id %in% shared_snps, ]

  # Run coloc
  eqtl_d <- list(
    beta = eqtl_region$beta, varbeta = eqtl_region$se^2,
    snp = eqtl_region$snp_id, position = eqtl_region$pos,
    type = "quant", MAF = eqtl_region$maf,
    N = max(eqtl_region$n), sdY = 1
  )
  gwas_d <- list(
    beta = gwas_region$beta, varbeta = gwas_region$se^2,
    snp = gwas_region$snp_id, position = gwas_region$pos,
    type = "quant", MAF = gwas_region$maf,
    N = max(gwas_region$n), sdY = 1
  )

  coloc.abf(dataset1 = eqtl_d, dataset2 = gwas_d)
}

# Example: iterate over GWAS lead SNPs
gwas_leads <- read.table("gwas_lead_snps.txt", header = TRUE)
coloc_results <- list()

for (i in seq_len(nrow(gwas_leads))) {
  res <- run_coloc_locus(
    gwas_leads$chr[i], gwas_leads$pos[i], gwas_leads$pos[i],
    eqtl_sumstats, gwas_sumstats
  )
  if (!is.null(res)) {
    coloc_results[[i]] <- res$summary
  }
}

# Filter for H4 > 0.75
coloc_h4 <- do.call(rbind, coloc_results)
coloc_positive <- coloc_h4[coloc_h4["PP.H4.abf", ] > 0.75, ]
cat("Loci with PP.H4 > 0.75:", nrow(coloc_positive), "\n")
```

## Sensitivity Analysis

```r
# Prior sensitivity: vary p12 to test robustness
priors <- c(1e-6, 1e-5, 5e-5, 1e-4)
sensitivity <- lapply(priors, function(p12) {
  coloc.abf(dataset1 = eqtl_data, dataset2 = gwas_data,
            p1 = 1e-4, p2 = 1e-4, p12 = p12)$summary
})

# Compare PP.H4 across priors
sapply(sensitivity, function(s) s["PP.H4.abf", ])
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| p1 (eQTL prior) | 1e-4 | ~1/10000 SNPs affect expression of a given gene |
| p2 (GWAS prior) | 1e-4 | ~1/10000 SNPs affect the trait |
| p12 (joint prior) | 1e-5 | p1 * p2 for independent traits; lower = more conservative |
| PP.H4 threshold | >= 0.75 | Standard for "strong evidence of shared causal variant" |
| PP.H4 suggestive | 0.5-0.75 | Worth follow-up validation |
| Min shared SNPs | 20 | Fewer SNPs = unreliable coloc |
| Window size | 500 kb | Include flanking LD region |

## Plant-Specific Notes

- **LD differences**: Plant LD decay varies greatly. Self-pollinating species
  (rice: ~100 kb, soybean: ~300 kb) have longer LD than outcrossers (maize: ~2 kb).
  Adjust coloc window accordingly. For long-LD species, consider conditional analysis
  to disentangle multiple independent signals.
- **GWAS from mixed models**: If GWAS used a mixed model (GAPIT, GEMMA), ensure
  effect sizes are on the same scale as the eQTL linear model. Standardize both
  phenotypes (sdY = 1).
- **Polyploid coloc**: Run coloc per subgenome. Cross-subgenome coloc is generally
  not meaningful due to homeolog expression complexity.
- **Multi-environment GWAS**: If GWAS comes from MET analysis, use per-environment
  summary statistics rather than BLUP GWAS. Environmental heterogeneity can create
  spurious coloc signals with tissue-specific eQTLs.
- **Low-accuracy genome**: If working with a scaffold-level assembly, coloc may be
  unreliable because gene and SNP positions on different scaffolds are not useful
  for colocalization.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "dataset1 and dataset2 must have same SNPs" | Non-overlapping SNP sets | Intersect by rsID or chr:pos before coloc |
| coloc.abf fails silently | Insufficient SNPs or monomorphic variants | Check MAF; ensure >= 20 shared SNPs |
| All loci PP.H4 < 0.5 | No true colocalization OR low power | Increase sample size; verify GWAS and eQTL hit same region |
| PP.H3 dominates (independent signals) | Different causal variants in LD | Perform conditional analysis; check LD pattern |
| coloc returns PP.H4 = 1.0 | Too few SNPs in region | Increase window; add more variants from imputation |
