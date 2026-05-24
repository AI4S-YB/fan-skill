# BLINK — Bayesian-information and Linkage-disequilibrium Iteratively Nested Keyway

**Goal:** High-power GWAS using LD-based marker binning with Bayesian model selection — fastest and most powerful method for large, high-density panels
**Approach:** LD-informed binning groups correlated SNPs → Bayesian information criterion (BIC) selects best SNP per bin → iterative fixed-effect testing with selected pseudo-QTNs as covariates
**Best for:** Large populations (>=500) with high-density markers (>=50K SNPs) and low-to-moderate stratification (lambda < 1.05)

## Prerequisites
- R 4.0+, GAPIT 3.0+
- Numeric genotype matrix (GAPIT numeric format: 0/1/2)
- High-density markers (>=50K recommended; BLINK's LD binning loses information at low density)
- Phenotype data matching genotype IDs

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(GAPIT3)

# ── Step 1: Load data ──
myGD <- read.table("${GENOTYPE_FILE}", header = TRUE, row.names = 1)
myGM <- read.table("${MARKER_MAP}", header = TRUE)
myY  <- read.table("${PHENOTYPE_FILE}", header = TRUE)

# ── Step 2: Quality control ──
# BLINK works best with high-quality, high-density data
maf_filtered  <- maf_filter(myGD, threshold = ${MAF_THRESHOLD})      # typically 0.05
miss_filtered <- missing_filter(myGD, threshold = ${MISS_THRESHOLD})  # typically 0.10

# ── Step 3: Run BLINK ──
myGAPIT <- GAPIT(
  G                  = myGD,
  GD                 = myGD,
  GM                 = myGM,
  Y                  = myY[, c(1, ${TRAIT_COL})],
  model              = "BLINK",
  PCA.total          = ${PCA_N},          # 3 if stratification is low, 5-10 if moderate
  Multiple_analysis  = FALSE,             # BLINK does NOT support Multiple_analysis
  SNP.fraction       = ${SNP_FRACTION},   # Fraction of SNPs to test, default 1 (all)
  kinship.cluster    = "average",         # Hierarchical clustering for kinship groups
  kinship.group      = "Mean",            # Mean kinship per group
  groupingFile       = NULL,              # Optional: pre-computed groups from FarmCPU
  threshold.output   = ${THRESHOLD_OUTPUT}  # Output threshold: 0.01 for top hits, 1 for all
)

# ── Step 4: Extract significant SNPs ──
# Output: GAPIT.BLINK.[trait].GWAS.Results.csv
results <- read.csv("GAPIT.BLINK.${TRAIT}.GWAS.Results.csv", header = TRUE)
sig_snps <- results[results$P.value < ${SIG_THRESHOLD}, ]  # Bonferroni: 0.05 / n_markers
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `PCA.total` | 3 | Moderate stratification (lambda 1.05-1.1) → 5; Admixed panels → 10 | BLINK's default PCA correction is less aggressive than FarmCPU's pseudo-QTN approach |
| `MAF_THRESHOLD` | 0.05 | Large panels (>=1000) → 0.01; Small panels (<300) → 0.05 | Larger panels can detect rarer variants reliably |
| `MISS_THRESHOLD` | 0.10 | Imputed data → 0.05; GBS data → 0.20 | Imputation quality affects the effective information content |
| `SNP_FRACTION` | 1 | Exploratory scan → 0.5 (faster); Full publication scan → 1 | Testing a subset speeds up initial runs but may miss signals |
| `THRESHOLD_OUTPUT` | 0.01 | Want all SNP results → 1; Only genome-wide significant → Bonferroni value | Controls output file size; 1 outputs every marker |
| `kinship.cluster` | "average" | Highly structured populations → "complete" | "complete" is more conservative for relatedness grouping |
| `SNP.test` | not set (all) | If `groupingFile` from prior FarmCPU run is available → set `SNP.test` count | Reduces computation by testing only grouped SNP representatives |

---

## When BLINK is Best

- **Large sample size (>=500)**: BLINK's Bayesian model selection benefits from more data
- **High marker density (>=50K SNPs)**: LD binning works better when bins are well-populated
- **Low population stratification (lambda < 1.05)**: BLINK's PCA correction is sufficient
- **You want maximum statistical power**: BLINK generally detects more true associations than FarmCPU or CMLM at equivalent FPR
- **Computational speed matters**: BLINK is the fastest GAPIT model for large datasets

### BLINK vs FarmCPU Decision Guide

| Criterion | Use BLINK | Use FarmCPU |
|-----------|:---:|:---:|
| Sample size | >= 500 | 200-500 |
| Marker density | >= 50K | 1K-50K |
| Stratification (lambda) | < 1.05 | > 1.05 |
| Speed | Fastest | Fast |
| False positive control | Good (BIC-based) | Very good (pseudo-QTN iteration) |
| Small datasets | Unstable | Works |
| Multiple_analysis support | No | Yes |

## When to Avoid

- **Low marker density (<10K SNPs)**: BLINK's LD-based binning loses too much information — fewer bins means lower resolution
- **Strong population stratification (lambda > 1.10)**: PCA correction alone is insufficient; use FarmCPU with its pseudo-QTN approach
- **Small samples (<200)**: Model selection becomes unstable; use CMLM or FarmCPU
- **Low-density legacy datasets** (SSR, RFLP): Use GLM or CMLM
- **Pedigree-based populations** (biparental, NAM with few families): Kinship-based mixed models are more appropriate

---

## Plant-Specific Notes

### Maize NAM (Nested Association Mapping) panels
- NAM populations have strong family structure — BLINK may need higher `PCA.total` (5-10)
- **Recommendation**: Run both BLINK and FarmCPU; compare overlap of significant SNPs
- Joint linkage-association mapping (JLAM) is preferred for NAM — BLINK alone may inflate false positives from family effects
- If using BLINK on NAM: first run FarmCPU with `Multiple_analysis = TRUE` to get `groupingFile`, then feed it to BLINK via the `groupingFile` parameter

### Large breeding populations (maize, rice, wheat diversity panels)
- BLINK excels here — large n + high density is the ideal use case
- For panels >2000 individuals: consider subsampling to validate that signals are not driven by subpopulation structure
- Compare BLINK results with FarmCPU on the full panel to identify consensus associations

### Polyploid species (wheat, cotton, canola)
- BLINK's LD binning may group homeologous SNPs from different subgenomes together if they share sequence similarity
- **Critical**: Ensure your marker map correctly assigns chromosomes/subgenomes to prevent cross-subgenome binning
- Run BLINK per subgenome as a validation step

### Self-pollinated crops with long-range LD (rice, soybean)
- Long LD blocks mean fewer independent bins — BLINK loses effective marker count
- FarmCPU may be more appropriate if LD extends >1 Mb

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Error: BLINK requires numeric genotype" | Genotype in HapMap/character format | Convert with `GAPIT.HapMap2numeric()` |
| Model runs but QQ plot shows severe deflation/inflation | `PCA.total` too low for the population structure | Increase `PCA.total` to 5-10; check lambda from GLM first |
| Manhattan plot has no peaks while FarmCPU shows peaks | BLINK's strict BIC penalizes large models — true associations may be penalized out | Use FarmCPU for final results; BLINK may be too conservative for this dataset |
| BLINK runs much slower than expected | `SNP_FRACTION` not set + very large dataset (>1M SNPs) | Set `SNP_FRACTION = 0.5` for initial scan; use LD-pruned subset |
| "kinship.group = Mean" error | All individuals assigned to single group | Check genotype quality; remove monomorphic markers; increase `kinship.cluster` granularity |
| Significant peaks at centromeric regions | LD binning groups non-causal centromeric SNPs with causal ones | Cross-reference with recombination rate map; flag centromeric peaks as potential false positives |
| Results differ substantially from FarmCPU | Different model assumptions — BLINK is more liberal, FarmCPU more conservative | Report both; treat overlapping peaks as high-confidence; divergent peaks require validation |
