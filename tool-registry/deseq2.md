# DESeq2 — Differential Expression Analysis

**Goal:** Differential expression analysis from count data using negative binomial GLM
**Approach:** Negative binomial GLM → Wald test → empirical Bayes shrinkage (apeglm)
**Best for:** Standard RNA-seq with ≥3 biological replicates per group

## Prerequisites
- R 4.0+, DESeq2 (Bioconductor), apeglm (for shrinkage)
- Raw count matrix (genes × samples) — NOT TPM/FPKM
- Sample metadata table (condition/group assignments)
- For transcript-level input: tximport + Salmon/Kallisto quantification

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(DESeq2)

# ── Step 1: Build DESeqDataSet ──
dds <- DESeqDataSetFromMatrix(
  countData = ${COUNTS},         # gene × sample matrix, raw counts
  colData   = ${METADATA},       # sample info, rownames must match colnames(counts)
  design    = ~ ${DESIGN}        # ~ condition (simple) / ~ batch + condition (batch-corrected)
)

# ── Step 2: Pre-filter low-count genes ──
# Keep genes with ≥ MIN_COUNT counts in at least MIN_SAMPLES samples
keep <- rowSums(counts(dds) >= ${MIN_COUNT}) >= ${MIN_SAMPLES}
dds <- dds[keep, ]

# ── Step 3: Set reference level ──
# Ensure the control group is the reference (first level)
dds$condition <- relevel(dds$condition, ref = "${CONTROL_GROUP}")

# ── Step 4: Run DESeq2 ──
dds <- DESeq(dds)

# ── Step 5: Extract results ──
res <- results(dds,
  contrast      = c("${VARIABLE}", "${TREAT}", "${CTRL}"),
  alpha         = ${ALPHA},           # FDR threshold, typically 0.05
  lfcThreshold  = ${LFC_THRESHOLD}    # 0 = no threshold; 1 = biological significance
)
res <- res[order(res$padj), ]
summary(res)

# ── Step 6: Shrink log2 fold changes (for visualization) ──
res_shrunk <- lfcShrink(dds, coef = "${COEF_NAME}", type = "apeglm")

# ── Step 7: Export ──
write.csv(as.data.frame(res), "${OUTPUT}_DESeq2_results.csv")
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `MIN_COUNT` | 10 | Low sequencing depth → 5 | Fewer reads per gene, lower threshold needed |
| `MIN_SAMPLES` | 3 (or n_reps) | **2 replicates → set to 2** | Otherwise all genes filtered out |
| `ALPHA` | 0.05 | Exploratory → 0.1; Strict publication → 0.01 | FDR stringency |
| `LFC_THRESHOLD` | 0 | **2 replicates → set to 1** | Reduces false positives from unstable variance estimates |
| `DESIGN` | ~ condition | With batch effects → ~ batch + condition | Avoids batch confounding |
| `ref` | control group | — | Always set explicitly, don't rely on alphabetical order |

---

## Replicate Count Strategy

| Reps | Method | Action |
|:---:|------|------|
| ≥ 3 | Standard DESeq2 | Optimal. Use default parameters. |
| 2 | DESeq2 with caution | **Must set `lfcThreshold=1`**. Power is low — expect fewer DEGs. Report honestly. |
| 1 | No statistical inference possible | Report fold change only. Do NOT report p-values. |

---

## Plant-Specific Notes

### Polyploid species (wheat, cotton, canola, potato)
- Homeologs have separate gene IDs (e.g., TraesCS1A01G0001, TraesCS1B01G0001, TraesCS1D01G0001)
- **Do NOT merge counts** — keep them as separate genes
- Bias analysis: compare DEG counts per subgenome after analysis

### Non-model species (apple, tea, most trees)
- No standard OrgDb for GO enrichment
- **Use Mercator (https://www.plabipd.de/portal/mercator4) for online GO annotation**, then import results
- Cross-reference with Arabidopsis homologs for functional inference
- If no GO results: it may be annotation poverty, not biological absence

### Tissue-specific zero-inflation
- Plant RNA-seq often has >50% zeros in tissue-specific genes — this is NORMAL
- DESeq2's independent filtering handles this automatically
- Do NOT filter genes that are expressed in ANY condition

### Circadian / time-of-day effects
- Sampling time matters. Circadian genes can vary >10-fold within 24h
- Note sampling time in metadata
- If time is uncontrolled, check whether top DEGs are known circadian genes

---

## Enrichment Strategy by Species

| Species category | GO enrichment approach |
|:---|------|
| Rice, maize, Arabidopsis | Plant Reactome first, then GO via OrgDb |
| Soybean, tomato, cotton | GO via OrgDb (if available) or AgriGO v2 |
| Apple, tea, most trees, non-model crops | **Mercator online annotation → GO** |
| Any species with >50 DEGs | If no enrichment found → report "may be annotation poverty" |

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "model matrix not full rank" | Too many covariates, not enough replicates | Simplify design formula; drop interaction terms |
| "all counts are zero" | Wrong gene annotation used for quantification | Check GFF/GTF file; gene IDs must match between annotation and count matrix |
| "contrast not found" | Condition name mismatch or wrong variable name | Check `levels(dds$condition)` or `colnames(colData(dds))` |
| Summary shows 0 DEGs | Power too low, or genuinely no differential expression | Check replicate count; try relaxing `lfcThreshold`; report honestly |
| padj = NA for all genes | Too few replicates to estimate dispersion | Minimum 2 replicates; if 2, use `lfcThreshold=1` |
| Error: "some values in assay are not integers" | Input is TPM/FPKM, not raw counts | Use `round()` to convert or re-quantify with featureCounts/HTSeq |
| rtracklayer::import() hangs on large GFF3 | GFF3 too large for single-threaded R parsing | Use Python `gzip + csv` or `gffutils` for fast extraction |
