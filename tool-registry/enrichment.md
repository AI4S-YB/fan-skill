# clusterProfiler — Functional Enrichment Analysis

**Goal:** Identify over-represented GO terms, KEGG pathways, and Plant Reactome pathways among DEGs, co-expression modules, or multi-omics features
**Approach:** Over-Representation Analysis (ORA) via hypergeometric test, with three strategies: (1) enrichGO with OrgDb for model species, (2) enricher with custom annotation for non-models, (3) Plant Reactome GMT for crop-specific pathways
**Best for:** After DE, WGCNA, GRN, or multi-omics integration — adds biological interpretation to gene lists

## Prerequisites
- R 4.0+, clusterProfiler (Bioconductor), org.*.eg.db (species-specific OrgDb)
- For non-model species: Mercator online annotation output, custom TERM2GENE mapping
- Gene list with species-appropriate gene IDs (RAP-DB, MSU, TAIR, or custom)

---

## Code Skeleton

This is NOT a fixed script. Choose the strategy matching your species category.

```r
library(clusterProfiler)
library(enrichplot)

# ── Step 1: Choose enrichment strategy based on species ──

# STRATEGY A: Model species with OrgDb (rice, Arabidopsis, maize, soybean, tomato)
enrich_result <- enrichGO(
  gene          = ${DEG_GENES},            # vector of gene IDs, e.g. c("LOC_Os01g01010", ...)
  OrgDb         = ${ORG_DB},               # org.Osativa.eg.db / org.At.tair.db / org.Sc.sgd.db
  keyType       = ${KEY_TYPE},             # "RAP-DB" → "RAPDB"; "MSU" → "MSU"; TAIR → "TAIR"; ENTREZ → "ENTREZID"
  ont           = "${ONT}",                # "BP" (Biological Process), "MF", "CC", or "ALL"
  pAdjustMethod = "${P_ADJUST}",          # "BH" (standard), "bonferroni" (strict)
  pvalueCutoff  = ${PVALUE_CUT},           # 0.05 (default), 0.01 (strict publication)
  qvalueCutoff  = ${QVALUE_CUT},           # 0.2 (exploratory), 0.05 (standard)
  minGSSize     = ${MIN_GS},               # 10 (standard), 5 (small gene lists)
  maxGSSize     = ${MAX_GS},               # 500 (standard)
  readable      = FALSE                    # TRUE to convert ENTREZID to symbol in output
)

# STRATEGY B: Non-model species — Mercator annotation → enricher
# Step B1: Upload genes to https://www.plabipd.de/portal/mercator4
# Step B2: Download Mercator GO annotation output (gene_id \t GO_term)
# Step B3: Build TERM2GENE mapping and run enricher

mercator_go <- read.table("${MERCATOR_OUTPUT}.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Expected columns: gene_id, GO_term (one row per gene-GO pair)
term2gene <- mercator_go[, c("GO_term", "gene_id")]

enrich_result <- enricher(
  gene          = ${DEG_GENES},
  TERM2GENE     = term2gene,
  pAdjustMethod = "${P_ADJUST}",
  pvalueCutoff  = ${PVALUE_CUT},
  qvalueCutoff  = ${QVALUE_CUT},
  minGSSize     = ${MIN_GS},
  maxGSSize     = ${MAX_GS}
)

# STRATEGY C: Plant Reactome ORA
# Download GMT for your species: https://plantreactome.gramene.org/download.html
gmt <- read.gmt("${PLANT_REACTOME}.gmt")   # e.g., "Osativa_plant_reactome.gmt"

enrich_result <- enricher(
  gene          = ${DEG_GENES},
  TERM2GENE     = gmt,
  pAdjustMethod = "${P_ADJUST}",
  pvalueCutoff  = ${PVALUE_CUT},
  qvalueCutoff  = ${QVALUE_CUT},
  minGSSize     = ${MIN_GS},
  maxGSSize     = ${MAX_GS}
)

# ── Step 2: KEGG enrichment (model species with KEGG annotation) ──
enrich_kegg <- enrichKEGG(
  gene          = ${ENTREZ_GENES},          # KEGG requires ENTREZ gene IDs
  organism      = "${KEGG_ORG}",             # "osa" (rice), "ath" (Arabidopsis), "zma" (maize)
  keyType       = "kegg",
  pAdjustMethod = "${P_ADJUST}",
  pvalueCutoff  = ${PVALUE_CUT},
  qvalueCutoff  = ${QVALUE_CUT}
)

# ── Step 3: Visualize results ──
dotplot(enrich_result, showCategory = 20, title = "${TITLE}")
barplot(enrich_result, showCategory = 20)
cnetplot(enrich_result, showCategory = 5, foldChange = ${FC_VECTOR})   # gene-concept network
heatplot(enrich_result, showCategory = 20, foldChange = ${FC_VECTOR}) # heatmap
emapplot(pairwise_termsim(enrich_result), showCategory = 30)          # enrichment map

# ── Step 4: Compare multiple gene lists ──
cmp <- compareCluster(
  geneClusters = list(
    up   = ${UP_GENES},
    down = ${DOWN_GENES}
  ),
  fun   = "enrichGO",
  OrgDb = ${ORG_DB},
  ont   = "${ONT}",
  pAdjustMethod = "${P_ADJUST}"
)
dotplot(cmp, showCategory = 15)

# ── Step 5: Export results ──
write.csv(as.data.frame(enrich_result), "${OUTPUT}_enrichment.csv")
```

---

## Species × Enrichment Method Decision Table

| Species category | Example species | Strategy | Tool / Pathway |
|:---|:---|:---|:---|
| **Model — high annotation** | Rice, Arabidopsis, maize | enrichGO → enrichKEGG | Plant Reactome GMT + GO (enrichGO) + KEGG (enrichKEGG) |
| **Model — moderate annotation** | Soybean, tomato, cotton | enrichGO + enricher | GO via OrgDb if available; Plant Reactome if GMT exists |
| **Non-model — draft annotation** | Apple, tea, peach, poplar | Mercator → enricher | Mercator for GO mapping; Plant Reactome GMT if available |
| **Non-model — no annotation** | Most trees, orphan crops | Mercator → enricher; AgriGO v2 fallback | Mercator online + AgriGO v2 (http://systemsbiology.cau.edu.cn/agriGOv2/) |
| **Any species with >50 DEGs** | — | If no enrichment found | Report "may be annotation poverty, not biological absence" |

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `pAdjustMethod` | "BH" | Bonferroni → "bonferroni" (stringent); exploratory → "none" (not recommended for publication) | BH controls FDR; Bonferroni is very conservative |
| `pvalueCutoff` | 0.05 | Strict publication → 0.01; Exploratory scan → 0.1 | Tighter threshold reduces false positives |
| `qvalueCutoff` | 0.05 | Few DEGs (<50) → 0.2; Large gene lists → 0.05 | Small gene lists produce fewer enriched terms; relax q-value to see signal |
| `minGSSize` | 10 | Small gene lists (<100) → 5; Very specific terms wanted → 3 | Too small → noisy GO terms; too large → misses specific biology |
| `maxGSSize` | 500 | Hyper-specific plant terms → 300 | Very large terms ("metabolic process") are usually uninformative |
| `ont` | "BP" | Mechanistic study → "BP" + "MF"; Localization study → "CC"; First pass → "ALL" | BP = regulation/process; MF = enzymatic function; CC = where |
| `keyType` | OrgDb-dependent | Rice: "RAPDB" vs "MSU" — pick the one matching your DEG IDs | ID mismatch = "No genes mapped" error |
| `readable` | FALSE | When output needs gene symbols for plotting | Converts ENTREZID to symbol; only works if OrgDb has symbol mapping |

---

## Non-Model Species Strategy

When enrichGO fails (no OrgDb exists for your species), follow this pipeline:

1. **Upload DEG list** to [Mercator4](https://www.plabipd.de/portal/mercator4) (accepts protein/nucleotide FASTA + gene ID list)
2. **Download GO annotation** from Mercator output (tab-separated: `gene_id \t GO_term`)
3. **Build TERM2GENE** mapping from Mercator output
4. **Run `enricher()`** with custom TERM2GENE (see Strategy B above)
5. **If no enrichment found**: report honestly — it may be annotation poverty, not biological absence. Cross-reference with Arabidopsis homologs for functional inference.
6. **Fallback**: [AgriGO v2](http://systemsbiology.cau.edu.cn/agriGOv2/) — web-based GO enrichment for 394 plant species, supports custom background

### Warning signs of annotation poverty (vs genuine absence of signal)
| Sign | Indicates |
|------|------|
| Corrected p > 0.05 for ALL GO terms | Could be annotation poverty |
| < 30% of input genes map to any GO term | Definitely annotation poverty — most genes lack annotation |
| Arabidopsis homolog of same genes shows strong enrichment | Confirms annotation poverty in target species |
| Enrichment present at p < 0.01 uncorrected but lost after correction | Likely genuine weak signal, not annotation poverty |

---

## Plant-Specific Notes

### Rice ID system — RAP-DB vs MSU (critical)
- Rice has **two incompatible gene ID systems**: RAP-DB (`Os01g0100100`) and MSU (`LOC_Os01g01010`)
- **enrichGO ID mapping will silently fail** if you provide RAP-DB IDs but keyType is set to MSU (or vice versa)
- Check your DEG IDs: RAP-DB = `OsXXg...` format; MSU = `LOC_OsXX...` format
- If your DE pipeline used RAP-DB: set `keyType = "RAPDB"` and use `org.Osativa.eg.db`
- If you need to convert: use `bitr()` with `fromType = "RAPDB"`, `toType = "MSU"`

### Plant Reactome GMT availability
- Download Plant Reactome GMT files: https://plantreactome.gramene.org/download.html
- Available species: Arabidopsis, rice, maize, soybean, tomato, grape, poplar, and more
- GMT format compatible with `enricher()` directly — no OrgDb needed
- Plant Reactome is **preferred over GO for crops** because pathways are curated specifically for plant biology
- Update frequency: check the Gramene download page for the latest release

### Gene family expansions in plants
- Plants have expanded gene families (e.g., NBS-LRR, cytochrome P450, MYB transcription factors)
- Enrichment of terms like "defense response" or "secondary metabolism" may reflect gene family architecture, not biology
- Use `simplify()` from clusterProfiler to collapse redundant GO terms (semantic similarity cutoff = 0.7)
- Cross-check enriched terms against genome-wide background of your species (some terms are always over-represented in plants)

### Polyploid species
- Homeologs in polyploids (wheat, cotton, canola) may have separate gene IDs per subgenome
- Do NOT merge homeologs for enrichment — keep them as separate input genes
- After enrichment, check whether enriched terms differ by subgenome (subgenome dominance pattern)
- If annotation exists for only one subgenome: report this as a limitation

### Time-of-day and circadian effects
- Circadian clock genes are highly enriched in GO categories related to "rhythm" and "photosynthesis"
- If your DEG list is unexpectedly enriched for these terms, verify that sampling time is consistent between conditions
- Arabidopsis circadian gene list available from: https://www.bioinf.manchester.ac.uk/circadian/

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No genes mapped" | Gene ID format mismatch (RAP-DB vs MSU for rice; TAIR10 vs TAIR9 for Arabidopsis) | Check ID format; use `bitr()` to convert between ID types; verify `keyType` matches your IDs |
| "No enrichment found" (all p > 0.05) | Too few DEGs, wrong background, or annotation poverty | Relax `qvalueCutoff` to 0.2; use all expressed genes as `universe`; if < 30% genes map → annotation poverty |
| "No enrichment found" (non-model species) | enrichGO used without OrgDb — fails silently | Switch to Mercator → enricher pipeline (Strategy B); verify OrgDb exists with `library(org.XX.eg.db)` |
| "universe is missing" or zero overlap | Background gene set doesn't contain any DEGs | Use all expressed genes (not whole genome) as `universe`; verify gene ID format matches between gene list and universe |
| "Error in bitr: fromType should be one of..." | Wrong keyType string for your species | Run `keytypes(org.Osativa.eg.db)` to list valid keyType values |
| "could not find function 'enrichGO'" | clusterProfiler not installed | `BiocManager::install("clusterProfiler")` |
| AnnotationDbi error (non-model species) | Trying to load OrgDb for species without one | Use Strategy B (Mercator + enricher); no OrgDb exists for non-model plants |
| "Error: TERM2GENE must be a data.frame" | GMT loaded incorrectly or Mercator output has wrong format | Verify `colnames(term2gene)` includes "term" and "gene" columns; use `read.gmt()` for Plant Reactome |
| enrichKEGG returns 0 results | KEGG organism code wrong or genes are not in KEGG | Verify organism code: `search_kegg_organism("${SPECIES}", by = "scientific_name")` |
| `simplify()` returns empty result | Semantic similarity threshold too strict | Raise cutoff from 0.7 to 0.8 or 0.9 |
| Duplicate enriched terms from gene family expansion | Redundant GO terms from large plant gene families | Use `simplify(enrich_result, cutoff = 0.7)` to collapse similar terms |
