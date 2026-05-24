# WGCNA — Weighted Gene Co-expression Network Analysis

**Goal:** Identify co-expressed gene modules, relate modules to phenotypic traits, and find hub genes driving module expression
**Approach:** Calculate gene-gene correlation matrix → raise to power beta (soft threshold) to achieve scale-free topology → hierarchical clustering → dynamic tree cut to define modules → correlate module eigengenes with traits
**Best for:** >=15 samples with continuous trait measurements; time-series or developmental gradient designs; tissue/organ atlases

## Prerequisites
- R 4.0+, WGCNA (CRAN), dynamicTreeCut
- Normalized gene expression matrix: genes as rows, samples as columns (WGCNA convention: columns = samples)
- Recommended input: variance-stabilized (VST) or rlog-transformed counts from DESeq2, or log2(TPM+1)
- Trait data table with matching sample identifiers (for module-trait association)
- **Minimum 15 samples** — fewer than 15 yields unreliable correlation estimates

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = ${N_THREADS})

# ── Step 1: Prepare expression data ──
# WGCNA expects samples as columns, genes as rows
datExpr <- t(${EXPR_MATRIX})  # Transpose if starting with genes-as-columns

# ── Step 2: Filter genes and samples ──
# Remove low-variance genes (they contribute noise)
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# Optional: keep top N most variable genes (for very large datasets >30K genes)
# var_genes <- apply(datExpr, 2, var)
# datExpr <- datExpr[, order(var_genes, decreasing = TRUE)[1:${TOP_N_GENES}]]

# ── Step 3: Choose soft-thresholding power ──
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "${NET_TYPE}",
                         verbose = 5)
# Plot to inspect:
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab = "Soft Threshold (power)", ylab = "Scale Free Topology R²")
# Choose power where R² > ${RSQ_CUTOFF} (recommended: 0.80-0.90)

# ── Step 4: Build network (blockwise for >5000 genes) ──
net <- blockwiseModules(
  datExpr        = datExpr,
  power          = ${POWER},            # From pickSoftThreshold
  networkType    = "${NET_TYPE}",       # "signed", "signed hybrid", or "unsigned"
  TOMType        = "${TOM_TYPE}",       # "signed" (recommended) or "unsigned"
  minModuleSize  = ${MIN_MODULE_SIZE},  # 30 for moderate, 50 for large datasets
  mergeCutHeight = ${MERGE_CUT_HEIGHT}, # 0.25 for moderate merging, 0.15 for conservative
  maxBlockSize   = ${MAX_BLOCK_SIZE},   # 5000-15000 depending on available RAM
  numericLabels  = FALSE,              # Use meaningful color names
  pamRespectsDendro = FALSE,
  saveTOMs       = FALSE,
  verbose        = 3
)

# ── Step 5: Module-trait association ──
# Correlate module eigengenes (MEs) with external traits
MEs <- net$MEs                       # Module eigengenes (samples × modules)
moduleTraitCor <- cor(MEs, ${TRAIT_DATA}, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# ── Step 6: Extract hub genes ──
# Genes with high module membership (kME) are hub candidates
kME <- signedKME(datExpr, MEs)
top_hubs <- apply(kME, 2, function(x) names(sort(x, decreasing = TRUE)[1:${N_HUBS}]))
moduleColors <- net$colors

# ── Step 7: Export results ──
write.csv(moduleTraitCor, "${OUTPUT}_moduleTrait_correlation.csv")
write.csv(cbind(gene = colnames(datExpr), module = moduleColors), "${OUTPUT}_gene_modules.csv")
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `POWER` | auto (pickSoftThreshold) | — | Always determine empirically; do NOT hard-code a power value |
| `RSQ_CUTOFF` | 0.80 | Small datasets (<20 samples) → 0.70; Strict publication → 0.90 | Lower n means noisier correlations — relaxed cutoff is pragmatic |
| `NET_TYPE` | "signed" | When direction of regulation is irrelevant → "unsigned" | Signed separates activation vs repression; unsigned groups both together — signed is typically more biologically meaningful |
| `MIN_MODULE_SIZE` | 30 | <5000 genes → 20; 5000-15000 genes → 30; >15000 genes → 50 | Scale with gene set size; too small = fragmented, too large = lumped |
| `MERGE_CUT_HEIGHT` | 0.25 | You want more distinct modules → 0.15; Want fewer, larger modules → 0.35 | Lower = more modules kept separate; affects biological interpretability |
| `MAX_BLOCK_SIZE` | 5000 | High-RAM server (>=64 GB) → 15000; Standard laptop (8-16 GB) → 3000 | Blockwise computes TOM on blocks; memory scales with block size squared |
| `deepSplit` | 2 (default) | More granular modules → 3-4; Fewer, broader modules → 0-1 | Controls sensitivity of dynamic tree cut; higher = more modules |

---

## When to Use blockwiseModules vs Standard WGCNA

| Criterion | blockwiseModules | Standard (manual block-by-block) |
|-----------|:---:|:---:|
| Number of genes | >5000 (recommended) | <5000 |
| RAM constraints | Scalable (block processing) | Entire TOM fits in memory |
| Ease of use | One function call | Manual control of each step |
| Customization | Moderate | Full control |
| Default choice | Yes for most cases | Only for very small datasets or manual debugging |

### Block size guidance by memory

| RAM available | maxBlockSize |
|:---|:---|
| 8 GB | 3000 |
| 16 GB | 5000 |
| 32 GB | 10000 |
| 64 GB | 20000 |
| >=128 GB | 30000 |

**Note**: If TOM computation exceeds memory, WGCNA will crash silently or the R session will be killed. Start conservatively with `maxBlockSize` and increase if resources allow.

---

## Plant-Specific Notes

### Tissue-specific co-expression modules
- Plant co-expression is strongly driven by tissue/organ identity — more so than in animals
- **Always include tissue/organ as a covariate** in trait association (regress out tissue effect before correlating module eigengenes with traits of interest)
- Tissue-specific modules (e.g., "seed-specific module," "root-specific module") are expected and biologically meaningful
- Modules expressed in multiple tissues often represent housekeeping or core metabolic functions — hub genes in these modules are candidates for pleiotropic regulators

### Developmental gradient modules
- Time-series experiments (seed development, fruit ripening, leaf senescence) produce gradient modules where eigengene expression correlates with developmental stage
- **Key**: check whether module eigengene smoothly increases/decreases with developmental time — these modules reflect developmental programs, not response to treatment
- For time-series: consider using `modulePreservation` across time points to identify transition points where network topology changes
- Ripening/fruit development: ethylene-responsive modules are nearly always present; verify by checking for ethylene-related TF hub genes (ERF, EIN3, etc.)

### Polyploid homeolog co-expression
- In allopolyploids (wheat, cotton, canola), homeologs (duplicated genes from different subgenomes) may:
  - **Co-express in the same module**: indicates conserved function / dosage balance
  - **Separate into different modules**: suggests subfunctionalization or neofunctionalization
- **Analysis tip**: Tag each gene with its subgenome of origin. After module assignment, test whether any module is enriched for a specific subgenome — subgenome-biased modules are evidence of subgenome dominance
- For autopolyploids (potato, alfalfa): homeologs are indistinguishable by sequence — module assignment may reflect allele-specific expression rather than homeolog divergence

### Stress response modules
- Abiotic stress (drought, salt, cold, heat) typically induces conserved modules enriched for:
  - Transcription factors: DREB, NAC, MYB, WRKY, bZIP families
  - Late embryogenesis abundant (LEA) proteins, heat shock proteins (HSPs)
  - ROS scavenging enzymes (SOD, CAT, APX)
- **Diagnostic check**: if >30% of DEGs from a DESeq2 contrast fall into 1-2 WGCNA modules, those modules are the transcriptional response to the treatment
- Biotic stress typically produces a distinct "defense module" enriched for NBS-LRR genes, PR proteins, and MAP kinase cascade components

### Circadian effects on co-expression
- Genes with strong circadian regulation can form artificial modules if samples are collected at different times of day
- If sampling time is uncontrolled: check whether top hub genes in each module are known circadian genes (CCA1, LHY, TOC1, GI, PRR family)
- For time-series across a day: consider removing circadian effect by including time-of-day as a covariate, or run WGCNA separately within each time point

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| No power achieves R² > 0.80 | Too few samples or noisy expression data | Accept R² > 0.70; increase sample size; filter low-variance genes more aggressively |
| "Error: cannot allocate vector of size X Mb" | `maxBlockSize` too large for available RAM | Reduce `maxBlockSize` (see Block size guidance table above) |
| All genes assigned to one "grey" (unassigned) module | `minModuleSize` too large or `deepSplit` too low | Reduce `minModuleSize` to 15-20; increase `deepSplit` to 3-4 |
| Too many modules (>50) with few genes each | `minModuleSize` too small or `mergeCutHeight` too low | Increase `minModuleSize` to 50; increase `mergeCutHeight` to 0.30-0.35 |
| Module eigengene has low correlation with trait (r < 0.3) | Module is not functionally related to the measured trait, or trait variation is too low | Report honestly; check whether any module correlates with an unexpected trait — this is often a discovery |
| Module eigengene correlates with batch/date instead of biological trait | Batch effects dominate expression variation | Include batch as covariate in trait table; or use `empiricalBayesLM` to regress out batch effects before WGCNA |
| Different modules show enrichment for the same GO terms | Modules represent sub-components of the same biological process; `mergeCutHeight` may be too low | Increase `mergeCutHeight` to 0.30; or report modules as sub-networks of the same process |
| Hub genes are mostly ribosomal proteins | Ribosomal protein genes are highly co-expressed and connected in nearly all tissues | Remove ribosomal protein genes before WGCNA if they dominate hub lists; they are a technical artifact for most biological questions |
| Signed network produces mostly negative eigengene-trait correlations | Signed network preserves direction — negative correlations are biologically meaningful and distinguish activated vs repressed modules | Do NOT use unsigned network to make them positive; interpret negative correlations as repression of that module |
| `goodSamplesGenes` removes all samples | Entire rows/columns are missing or constant | Check input data: missing values must be removed, constant-expression genes must be filtered |
| Module preservation test shows all modules are not preserved in test dataset | Test dataset is biologically too different (different tissue, species, or treatment) | Verify that comparison is biologically meaningful; some loss of preservation is expected across tissues |
