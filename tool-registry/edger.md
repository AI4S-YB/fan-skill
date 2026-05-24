# edgeR — Differential Expression Analysis

**Goal:** Differential expression analysis from count data using empirical Bayes quasi-likelihood F-tests
**Approach:** DGEList → filterByExpr → TMM normalization → estimateDisp → glmQLFit → glmQLFTest
**Best for:** Low-replicate designs (2 replicates), or when TMM normalization is preferred for plant RNA composition bias

## Prerequisites
- R 4.0+, edgeR (Bioconductor)
- Raw count matrix (genes x samples) — NOT TPM/FPKM
- Sample metadata with group/condition assignments
- For transcript-level input: tximport + Salmon/Kallisto quantification

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(edgeR)

# ── Step 1: Create DGEList ──
y <- DGEList(
  counts = ${COUNTS},           # gene × sample matrix, raw integer counts
  group  = ${GROUP}             # factor defining sample groups (e.g., group <- factor(c("WT","WT","mut","mut")))
)

# ── Step 2: Filter low-expression genes ──
# filterByExpr keeps genes with sufficient counts for statistical analysis
keep <- filterByExpr(
  y,
  min.count      = ${MIN_COUNT},      # default 10; lower for shallow sequencing
  min.total.count = ${MIN_TOTAL}      # default 15
)
y <- y[keep, , keep.lib.sizes = FALSE]

# ── Step 3: Normalize (TMM) ──
y <- calcNormFactors(y, method = "TMM")

# ── Step 4: Build design matrix ──
design <- model.matrix(~ ${DESIGN})   # ~ 0 + group for explicit contrasts, or ~ group for treatment-vs-control

# ── Step 5: Estimate dispersion ──
# For 2+ replicates: use estimateDisp (QL pipeline)
# For 0 replicates: see "No Replicates" section below
y <- estimateDisp(y, design)

# ── Step 6: Fit quasi-likelihood model ──
fit <- glmQLFit(y, design, robust = TRUE)  # robust = TRUE helps with outlier genes

# ── Step 7: Test ──
qlf <- glmQLFTest(fit, contrast = ${CONTRAST})  # e.g., c(-1, 1) for two-group comparison
top <- topTags(qlf, n = Inf, sort.by = "PValue")

# ── Step 8: Export ──
write.csv(as.data.frame(top), "${OUTPUT}_edgeR_results.csv")
summary(decideTests(qlf, p.value = ${ALPHA}))
```

### No Replicates

```r
# When you have NO replicates — USE WITH EXTREME CAUTION
# Cannot estimate biological variability. Report fold change only.
bcv <- 0.4  # Assumed biological CV; 0.4 for well-controlled organisms, 0.1 for cell lines
et <- exactTest(y, dispersion = bcv^2)
top <- topTags(et, n = Inf)
# Do NOT report p-values in publications without replicates
```

### Paired / Block Design

```r
# When samples are paired (e.g., same plant before/after treatment)
design <- model.matrix(~ ${BLOCK} + ${TREATMENT})
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = ${TREATMENT_COEF})
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `MIN_COUNT` | 10 | Low sequencing depth → 5 | Fewer reads per gene, lower threshold needed |
| `method` | "TMM" | Extreme library composition bias → "RLE" or "upperquartile" | TMM is robust but can fail with dominant RNA species |
| `robust` | TRUE | With many outliers → TRUE; Small n → FALSE | Robust estimation helps but can be unstable with few samples |
| `design` | ~ group | Paired design → ~ block + group; Multi-factor → ~ batch + treatment | Controls for nuisance variation |
| `contrast` | c(-1, 1) or coef | Any pairwise comparison via contrast vector | Must match design matrix column order |
| `ALPHA` (FDR) | 0.05 | Exploratory → 0.1; Strict publication → 0.01 | FDR stringency |
| `bcv` (no reps) | 0.4 | Plants with high genetic variation → 0.6; Clonal plants → 0.2 | Higher CV for genetically diverse populations |

---

## Replicate Count Strategy

| Reps | Method | Action |
|:---:|------|------|
| ≥ 3 | Standard edgeR QLF | Optimal. Use `robust = TRUE` in glmQLFit. |
| 2 | edgeR QLF with caution | edgeR excels here. QL F-test outperforms Wald test at low n. Still: power is low, expect fewer DEGs. Report honestly. |
| 1 | No statistical inference possible | Use `exactTest()` with assumed BCV. Report fold change only. Do NOT report p-values. |
| 0 technical | No actual replicates | Pooled technical reps do NOT count as biological replicates. Treat as n=1. |

---

## When to Choose edgeR over DESeq2

| Scenario | Recommendation | Rationale |
|----------|:---:|------|
| 2 biological replicates | **edgeR** | QL F-test has better type I error control with small n than Wald test |
| Strong RNA composition bias (seed vs leaf) | **edgeR** | TMM normalization explicitly handles compositional differences |
| No replicates (exploratory) | **edgeR** (exactTest) | DESeq2 requires at least 2 replicates to run; edgeR can proceed with assumed BCV |
| ≥ 3 replicates, standard design | DESeq2 | Slightly more conservative; out-of-memory prediction handles low counts better |
| Multi-factor design with interactions | DESeq2 | Easier formula interface; edgeR requires manual design matrix construction |
| Transcript-level input (Salmon/Kallisto) | Either | Both integrate with tximport; bias correction available in both |

---

## Plant-Specific Notes

### TMM for seed vs leaf (very different RNA composition)
- Seed transcriptomes are dominated by storage protein RNAs (up to 80% of mRNA)
- Leaf transcriptomes have high Rubisco (rbcL/rbcS) mRNA
- **TMM explicitly models the "most genes are not differentially expressed" assumption** — this holds even with extreme compositional differences
- If TMM fails (M-values show bimodal distribution), try `method = "upperquartile"` or filter dominant transcripts (>50% of library) before normalization
- DESeq2's median-of-ratios can be distorted by dominant transcripts; edgeR's trimmed mean is more robust

### Polyploid species (wheat, cotton, canola, potato)
- Homeologs have separate gene IDs — **do NOT merge counts across subgenomes**
- Keep each homeolog as a separate row in the count matrix
- TMM normalization is generally safe: it trims the most extreme log-ratios (top/bottom 30% by default), so polyploid-specific expression differences are accommodated
- After analysis: compare DEG counts per subgenome to check for homeolog expression bias
- For allopolyploids (e.g., wheat A/B/D, cotton At/Dt), the "different RNA composition" issue is even more pronounced — TMM is well-suited

### Tissue-specific zero-inflation
- Plant RNA-seq often has >50% zeros in tissue-specific genes — this is NORMAL
- `filterByExpr()` handles this automatically by requiring expression in a minimum number of samples
- Do NOT manually filter genes that are expressed in ANY condition

### Circadian / time-of-day effects
- Sampling time matters. Circadian genes can vary >10-fold within 24h
- Note sampling time in metadata
- If time is uncontrolled, check whether top DEGs are known circadian genes

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "no residual df" | No replicates remaining after design | Minimum 2 biological replicates needed. If genuinely n=1, use `exactTest()` with assumed BCV — report fold change only |
| "contrast vector wrong length" | Length mismatch between contrast and design matrix columns | Check `colnames(design)`; contrast length must equal number of coefficients |
| "NA/NaN in results" | Genes with zero counts in all samples of one group | These genes are filtered by `filterByExpr()` — if many remain, lower `min.count` |
| "Dispersion estimate failed" | Too few replicates or extreme outliers | Try `estimateGLMCommonDisp()` first, then `estimateGLMTagwiseDisp()`; or remove outlier samples |
| filterByExpr filters ALL genes | `min.count` too high or all samples have very low depth | Lower `min.count` to 5; check library sizes with `y$samples$lib.size` |
| BCV estimate wildly different from DESeq2 | Different dispersion estimation methods | Normal: edgeR estimates genewise, DESeq2 shrinks toward trend. Reconcile by checking `plotBCV(y)` |
| CPM vs count confusion | `DGEList()` received CPM/TPM instead of raw counts | edgeR requires raw integer counts. If only normalized data is available, use `limma-voom` instead |
