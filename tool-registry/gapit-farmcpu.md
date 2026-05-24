# GAPIT FarmCPU — GWAS with Fixed and Random Model Circulating Probability Unification

**Goal:** GWAS with iterative fixed-effect testing + random-effect polygenic background — high statistical power with low false positives
**Approach:** Alternating fixed-effect model (test markers one-by-one) → random-effect model (estimate polygenic background from associated markers) → iterate until no new associations
**Best for:** Medium-to-large populations with population stratification (lambda > 1.05) and marker density >= 1000 SNPs

## Prerequisites
- R 4.0+, GAPIT 3.0+
- Numeric genotype matrix (GAPIT numeric format: 0/1/2 or 0/1 for heterozygous)
- Phenotype data with genotype IDs matching the genotype matrix
- Marker density >= 1000 SNPs (FarmCPU needs enough markers to model polygenic background)
- Kinship matrix optional — FarmCPU estimates its own pseudo-kinship internally

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(GAPIT3)

# ── Step 1: Load and QC genotype data ──
# GAPIT numeric format: rows = taxa, columns = markers, values = 0/1/2
myGD <- read.table("${GENOTYPE_FILE}", header = TRUE, row.names = 1)
myGM <- read.table("${MARKER_MAP}", header = TRUE)   # SNP, Chromosome, Position
myY  <- read.table("${PHENOTYPE_FILE}", header = TRUE)

# Filter by MAF and missing rate (adjust thresholds to your data)
maf_filtered  <- maf_filter(myGD, threshold = ${MAF_THRESHOLD})  # typically 0.05
miss_filtered <- missing_filter(myGD, threshold = ${MISS_THRESHOLD})  # typically 0.10

# ── Step 2: Run FarmCPU ──
myGAPIT <- GAPIT(
  G                  = myGD,              # Genotype matrix (numeric)
  GD                 = myGD,              # Same as G (GAPIT convention)
  GM                 = myGM,              # Marker map
  Y                  = myY[, c(1, ${TRAIT_COL})],  # Taxon + trait column(s)
  model              = "FarmCPU",
  PCA.total          = ${PCA_N},          # 3-5; controls population structure
  Multiple_analysis  = TRUE,              # Multiple iterations improve stability
  model.selection    = TRUE,              # Report pseudo-QTN count
  bin.size           = ${BIN_SIZE},       # Window size in bp for binning, default 1e6
  bin.selection      = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
)

# ── Step 3: Interpret outputs ──
# Key files written to working directory:
#   GAPIT.FarmCPU.[trait].GWAS.Results.csv  → SNP, Chr, Pos, P.value, maf, effect
#   GAPIT.FarmCPU.[trait].Manhattan.png      → Manhattan plot
#   GAPIT.FarmCPU.[trait].QQ-Plot.png        → QQ plot
#   GAPIT.FarmCPU.Multiple_analysis.pdf      → Convergence across iterations
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `PCA.total` | 3 | Strong stratification (lambda > 1.2) → 5; Admixed panels → 5-10 | More PCs capture finer population structure |
| `MAF_THRESHOLD` | 0.05 | Small populations (<100) → 0.02; High-density arrays → 0.01 | Fewer individuals need lower MAF to retain enough markers |
| `MISS_THRESHOLD` | 0.10 | Low-coverage GBS → 0.20 | GBS data inherently has more missing calls |
| `BIN_SIZE` | 1e6 (1 Mb) | Small genome (<500 Mb) → 5e5; Large genome (>3 Gb) → 5e6 | Bin size should capture LD decay distance for the species |
| `Multiple_analysis` | TRUE | Exploratory (quick look) → FALSE | Without it, convergence cannot be assessed |
| `model.selection` | TRUE | Always keep TRUE | Shows how many pseudo-QTNs were used — diagnostic for model fit |
| `bin.selection` | default sequence | Small genome → use shorter sequence starting at 5 | Controls the range of bin sizes tested for optimal power |

---

## When FarmCPU Shines

- **Population stratification present (lambda > 1.05)** but you want higher power than CMLM
- **Marker density adequate (>1000 SNPs)** — FarmCPU needs markers to estimate polygenic effects
- **Moderate to large sample size (>=200)** — iterative model converges reliably
- **You want fewer false positives than GLM** but more power than strict MLM/CMLM

### FarmCPU vs CMLM Decision Guide

| Criterion | Use FarmCPU | Use CMLM |
|-----------|:---:|:---:|
| Population structure (lambda) | > 1.05 | > 1.05 (both work) |
| Marker density | >= 1000 SNPs | Any |
| Sample size | >= 200 | >= 100 |
| Computational speed | Faster | Slower |
| False positive control | Good | Conservative |
| Statistical power | High | Moderate |

## When to Avoid

- **Low marker density (<1000 SNPs)**: FarmCPU cannot properly model the polygenic background — use CMLM or GLM instead
- **Very small samples (<100)**: Model fitting becomes unstable; pseudo-QTN selection is unreliable
- **Low-density legacy datasets** (e.g., SSR markers, RFLP): Switch to CMLM with kinship matrix
- **Strong family structure with few families**: Mixed model with kinship performs better

---

## Plant-Specific Notes

### Polyploid species (wheat, cotton, canola, potato, sugarcane)
- FarmCPU handles polyploids better than most GWAS methods because pseudo-QTNs absorb background signals from all subgenomes simultaneously
- **Subgenome-specific analysis**: Run FarmCPU separately per subgenome to check for subgenome-biased associations — but report the whole-genome result as primary
- **Dosage scoring matters**: For allopolyploids, ensure genotypes are scored correctly (0/1/2 dosage, not presence/absence)
- **Stability check**: Run with `Multiple_analysis = TRUE` and verify that significant peaks are consistent across iterations — polyploid peaks can shift slightly

### Self-pollinated crops (rice, soybean, wheat, tomato)
- Population structure tends to be stronger — increase `PCA.total` to 5
- LD decay is typically longer — increase `BIN_SIZE` to 2-5 Mb

### Cross-pollinated crops (maize, sunflower, rye)
- LD decays faster — use default or smaller `BIN_SIZE` (5e5)
- Generally less stratification — `PCA.total = 3` is sufficient

### Perennial crops (apple, grape, tea, poplar)
- Often have high heterozygosity — ensure GAPIT numeric coding handles heterozygotes correctly
- Clonal propagation = identical genotypes represented by a single individual in the panel
- Use a kinship matrix from pedigree or markers if clonal relationships are complex

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Model did not converge" | Insufficient markers (<1000) or samples (<100) for pseudo-QTN selection | Switch to CMLM or BLINK; reduce `PCA.total` |
| "Too few markers after filtering" | MAF or missing rate thresholds too strict | Relax `MAF_THRESHOLD` to 0.02 or `MISS_THRESHOLD` to 0.20 |
| Manhattan plot shows no peaks above threshold | True absence of large-effect loci, or power too low | Check sample size; try Bonferroni threshold vs FDR; report negative result honestly |
| QQ plot shows severe inflation (lambda >> 1.1) | Residual population structure or cryptic relatedness | Increase `PCA.total`; add kinship matrix; check for outlier individuals via PCA |
| "Error in GAPIT: genotype matrix not numeric" | Genotype data in HapMap or VCF format | Convert to GAPIT numeric format (0/1/2) using `GAPIT.HapMap2numeric()` or `GAPIT.Numerical2numeric()` |
| Multiple_analysis plot shows unstable peaks | Model oscillating between different pseudo-QTN sets | Increase `bin.selection` range; verify data quality (no sample swaps) |
| Significant peaks disappear after adding PCs | Peaks were false positives from population structure | These were not true associations — report the corrected result |
