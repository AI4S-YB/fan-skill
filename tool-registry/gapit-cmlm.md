# GAPIT CMLM (Compressed Mixed Linear Model)

**Goal:** GWAS with compressed MLM — controlling both population structure and kinship simultaneously
**Approach:** Import genotypes → PCA for population structure → compute kinship (K matrix) → compress groups → test each SNP with MLM → Manhattan/Q-Q output
**Best for:** Inbred species, small populations (<200), low-to-medium density markers (<50K SNPs)

## Prerequisites
- R 4.0+, GAPIT 3.0+
- Genotype: Hapmap (preferred) or numerical format
- Phenotype: CSV/TXT, first column = sample ID (Taxa), subsequent columns = trait values
- Sample IDs must match exactly between genotype and phenotype files

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(GAPIT3)

# ── Step 1: Load data ──
myY  <- read.csv("${PHENOTYPE_FILE}", head = TRUE)   # Taxa, trait1, trait2, ...
myG  <- read.delim("${GENOTYPE_HAPMAP}", head = FALSE) # Hapmap format
myGM <- myG[, c(1, 3, 4)]                              # SNP name, chromosome, position

# ── Step 2: Run CMLM ──
myGAPIT <- GAPIT(
  G                = myG,              # Genotype (Hapmap)
  Y                = myY[, 1:2],       # Phenotype: Taxa + ONE trait at a time
  model            = "CMLM",

  # ── Population structure (PCA) ──
  PCA.total        = ${PCA_TOTAL},     # 3-5 for most plant populations
  PCA.col          = NULL,             # Use first N PCs; NULL = auto
  PCA.3d           = FALSE,            # 3D PCA plot (set TRUE for publication figures)

  # ── Kinship / compression ──
  kinship.cluster  = "${KINSHIP_CLUSTER}",  # "average", "complete", "ward"
  kinship.group    = "Mean",                # Kinship type: "Mean" (default) or "Max"
  group.from       = ${GROUP_FROM},         # Start of compression grid
  group.to         = ${GROUP_TO},           # End of compression grid
  group.by         = ${GROUP_BY},           # Step size (e.g., 100)

  # ── Output control ──
  cutOff           = ${CUTOFF},         # P-value threshold for Manhattan output (0.05 default)
  file.output      = TRUE               # Set FALSE for interactive exploration
)

# ── Step 3: Extract significant SNPs ──
# Results are in myGAPIT$GWAS — a list with one element per trait
sig_snps <- myGAPIT$GWAS[[1]][myGAPIT$GWAS[[1]]$P.value < ${SIG_THRESHOLD}, ]
write.csv(sig_snps, "${OUTPUT}_significant_SNPs.csv", row.names = FALSE)
```

### Numeric Genotype Format

```r
# If data is in numerical format (0/1/2 coding):
myGD <- read.csv("${NUMERIC_GENOTYPE}", head = TRUE)  # Taxa + SNP columns
myGM <- read.csv("${SNP_MAP}", head = TRUE)           # SNP, Chr, Pos
myGAPIT <- GAPIT(
  GD = myGD,
  GM = myGM,
  Y  = myY[, 1:2],
  model = "CMLM",
  PCA.total = ${PCA_TOTAL},
  kinship.cluster = "${KINSHIP_CLUSTER}",
  group.from = ${GROUP_FROM},
  group.to = ${GROUP_TO},
  group.by = ${GROUP_BY}
)
```

### Multi-Trait Run

```r
# GAPIT supports running multiple traits in one call:
myGAPIT <- GAPIT(
  G = myG,
  Y = myY,           # Include ALL trait columns
  model = "CMLM",
  PCA.total = ${PCA_TOTAL}
)
# myGAPIT$GWAS is a list with results for each trait
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `PCA.total` | 3 | Strong population structure (λ > 1.2) → 5; Admixed population → 5-10; Low structure → 0-2 | Balances confounding control vs over-correction |
| `kinship.cluster` | "average" | Many small clusters expected → "ward"; Outliers present → "complete" | Average linkage is standard; Ward gives more balanced groups |
| `group.from` | max(50, n/4) | Very small panel (<80) → n/2; Large panel (>500) → 100 | Start compression from a reasonable fraction of population |
| `group.to` | n (all samples) | For speed with n > 1000 → min(n, 1000) | Full compression is ideal but computationally expensive |
| `group.by` | max(10, n/20) | n < 100 → 5; n > 500 → 50 | Finer grid for small populations, coarser for large |
| `cutOff` | 0.05 | Exploratory → 0.01; Full table → 1.0 | Lower cutoff filters output files; use 1.0 to get all results |
| `file.output` | TRUE | Interactive exploration → FALSE | Suppresses file writing for rapid iteration |

---

## When to Choose CMLM vs FarmCPU vs BLINK

| Scenario | Recommended Model | Rationale |
|----------|:---:|------|
| **Small population (<200)** | **CMLM** | FarmCPU/BLINK need many markers to estimate polygenic background reliably |
| **Low-density markers (<50K)** | **CMLM** | CMLM tests each marker individually; FarmCPU/BLINK bin markers and need density |
| **Inbred species** (rice, wheat, soybean, cotton) | **CMLM** | K matrix captures relatedness well in structured inbred panels |
| **High genomic inflation (λ > 1.2)** | FarmCPU | FarmCPU's iterative approach explicitly separates false positives from population structure |
| **Large population (≥500)** | FarmCPU or BLINK | These methods scale better and have higher power in large panels |
| **High-density markers (≥50K)** | BLINK | BLINK's LD-based binning requires high density to define bins |
| **Moderate population (200-500), moderate density** | FarmCPU | Good balance of power and false positive control |
| **Outcross species** (maize, sunflower, perennial ryegrass) | FarmCPU or BLINK | LD decays faster in outcrossers, favoring binning approaches |
| **First pass / exploratory** | CMLM | Simpler output, easier to diagnose problems |
| **Replication / validation** | Match original study | Consistency with prior work aids interpretation |

### Cross-Reference
- **FarmCPU**: See `tool-registry/gapit-farmcpu.md` — iterative fixed + random model; better for structured/admixed populations
- **BLINK**: See `tool-registry/blink.md` — LD-binned GWAS; best for large high-density panels

---

## Replicate / Population Size Strategy

| Population size | Recommendation | Action |
|:---:|------|------|
| < 100 | CMLM only | Power is inherently limited. Focus on large-effect loci. Report Q-Q plot honestly. |
| 100-200 | CMLM (primary) | Optimal range for CMLM. Compression works well. |
| 200-500 | CMLM or FarmCPU | Run CMLM first. If λ > 1.1, also run FarmCPU and compare. |
| 500-2000 | FarmCPU or BLINK | CMLM becomes computationally slow and less powerful. |
| > 2000 | BLINK | Large populations with dense markers favor BLINK. |

---

## Plant-Specific Notes

### Inbred species (rice, wheat, soybean, cotton, barley, most self-pollinating crops)
- CMLM is the **default recommendation** for initial GWAS
- K matrix is essential — inbred panels often have cryptic relatedness that PCA alone cannot capture
- **Hapmap format**: Numeric genotype coding (0/1/2 for major allele count at bi-allelic SNPs)
- For inbred species, heterozygosity should be near zero — check Hapmap quality first
- Inbred panels often show strong population structure; always report the Q-Q plot and genomic inflation factor (λ)

### Outcross species (maize, sunflower, oil palm, perennial ryegrass, most trees)
- CMLM can still work but **PCA.total should be higher (5-10)** to control population structure
- Consider **FarmCPU or BLINK** if the population is large and markers are dense
- Outcrossers have more heterozygosity — ensure Hapmap encoding handles heterozygous calls correctly
- LD decays much faster in outcrossers, so CMLM's single-marker approach may miss small-effect QTL in bins
- Action: Run CMLM first. If λ > 1.15 or Manhattan shows severe inflation, switch to FarmCPU.

### Polyploid GWAS warning
- **CMLM (and GAPIT in general) assumes diploid genotype coding**. Using it on polyploid species without proper allele dosage models is problematic:
  - **Allotetraploids** (wheat, cotton, canola): Can use diploid models IF subgenome-specific markers are used and analyzed separately per subgenome
  - **Autopolyploids** (potato, alfalfa): Standard GWAS tools (including GAPIT) are generally NOT appropriate. Use tools like GWASpoly, SHEsis, or polyRAD that model allele dosage
- If you must use GAPIT on a polyploid: analyze each subgenome separately with subgenome-specific markers; report this limitation clearly

### Centromere false positives
- Plant centromeric regions often show **spurious GWAS peaks** due to suppressed recombination and high LD, NOT true associations
- After running CMLM, check whether top SNPs map to known centromeric/pericentromeric regions
- For species with annotated centromeres (rice, maize, Arabidopsis, soybean): filter or flag centromeric SNPs
- For non-model species: flag clusters of extremely significant SNPs spanning large physical distances as potential centromeric artifacts

### Low-density markers (GBS, RAD-seq)
- CMLM handles low-density markers **better than FarmCPU/BLINK** because it does not require LD-based binning
- With GBS data (5K-20K markers), CMLM is the safe default
- Be aware that low-density GWAS will miss small-effect QTL and may only detect large-effect loci
- Missing data rate matters: impute or filter SNPs with >20% missing calls before running GAPIT

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "GAPIT not found" | Package not installed | `devtools::install_github("jiabowang/GAPIT3")` — requires `devtools` package |
| "subscript out of bounds" | Mismatch between genotype and phenotype sample IDs | Check `intersect(myG[,1], myY[,1])`; ensure Taxa IDs match exactly (case, whitespace, special characters) |
| "NA/NaN in phenotype" | Missing values in phenotype column | Remove rows with `myY <- myY[!is.na(myY[,2]), ]` or impute; GAPIT cannot handle missing phenotypes |
| "Error in eigen" | Singular kinship matrix — identical or highly related individuals | Remove duplicate or clonal samples; check for identical genotypes |
| Manhattan plot is flat (no peaks) | Genuinely no significant associations, OR too stringent threshold | Check Q-Q plot first; if λ ≈ 1 and no deviation, trait may genuinely lack large-effect loci |
| Q-Q plot shows early deviation but flat Manhattan | Polygenic trait with many small-effect loci | Expected for complex traits. Report honestly. Consider multi-locus methods. |
| "cannot open file" for output | File permissions or path issues | Set `file.output = FALSE` for interactive work; check write permissions in working directory |
| Hapmap header mismatch | Row 1 of Hapmap expected to contain "rs#" header | Verify Hapmap format: row 1 = "rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID", "QCcode", then sample IDs |
| GAPIT hangs on large datasets | K matrix computation is O(n²) in memory | Reduce sample size or use `kinship.algorithm = "EMMA"` for approximation; consider FarmCPU for large panels |
| Centromeric peak with high significance | Suppressed recombination creates LD blocks spanning Mb | Cross-reference top SNPs with known centromere positions; if unavailable, flag clusters with extremely high -log10(p) spanning >5 Mb |
