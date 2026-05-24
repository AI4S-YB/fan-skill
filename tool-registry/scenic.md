# SCENIC — Single-Cell Regulatory Network Inference and Clustering

**Goal:** From expression matrix to regulons — infers co-expression modules between
TFs and targets (via GENIE3/GRNBoost2), prunes indirect targets via cis-regulatory
motif enrichment (RcisTarget), and scores regulon activity per sample (AUCell).

**Best for:** >= 20 samples, TF annotation available, regulon activity needed.

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| minGenes | 20 | Small TF families (<10 members): lower to 10; large genomes (wheat, maize): raise to 50 | Small modules with few genes are less likely to represent real regulons; large genomes produce more spurious co-expression |
| motif database | Species-specific (cisTarget) | Non-model plant: use Arabidopsis motifs as proxy; well-annotated model (Arabidopsis, rice): use DAP-seq + ChIP-seq combined | Plant motif databases are incomplete; Arabidopsis motifs can identify conserved binding preferences across angiosperms |
| nes_threshold (ctx) | 3.0 | Low recovery (<50 regulons): reduce to 2.5; excess regulons (>500): raise to 4.0 | Normalized enrichment score controls the stringency of motif-based pruning; non-model species often need relaxed thresholds |
| auc_threshold (binarization) | 0.05 | Sparse data (<100 cells): use 0.10; dense scRNA-seq: auto-detect via bimodal distribution | Binarization threshold depends on signal-to-noise ratio; sparse data has lower AUC values overall |

## Prerequisites

- R 4.0+ or Python 3.8+ (pySCENIC)
- R packages: SCENIC, RcisTarget, AUCell, GENIE3 (Bioconductor)
- Python packages: pyscenic, arboreto, ctx-core (cisTarget databases)
- cisTarget motif databases (download from https://resources.aertslab.org/cistarget/)
- Species-specific motif rankings (.feather) and motif-TF annotation (.tbl)

## Basic Usage (pySCENIC)

```bash
# Step 1: GRN inference (GRNBoost2)
pyscenic grn \
    expression_matrix.tsv \
    tf_genes.txt \
    -o grn_adjacencies.tsv \
    --num_workers 4

# Step 2: Regulon pruning (cisTarget)
pyscenic ctx \
    grn_adjacencies.tsv \
    motifs-rankings.feather \
    --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname expression_matrix.tsv \
    --output regulons.tsv \
    --num_workers 4

# Step 3: AUCell — regulon activity per sample
pyscenic aucell \
    expression_matrix.tsv \
    regulons.tsv \
    -o aucell_output.tsv \
    --num_workers 4
```

## Basic Usage (R SCENIC)

```r
library(SCENIC)

# Initialize
scenicOptions <- initializeScenic(
    org = "plant",                     # Use plant-specific settings
    dbDir = "cisTarget_databases",
    nCores = 4
)

# Load expression
exprMat <- as.matrix(read.table("expression_matrix.tsv",
    header = TRUE, row.names = 1, sep = "\t"))
genesKept <- geneFiltering(exprMat, scenicOptions)

# Run GENIE3 -> RcisTarget -> AUCell
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)
runSCENIC_4_aucell_binarize(scenicOptions)

# Regulon specificity score (RSS)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC = getAUC(regulonAUC),
               cellAnnotation = cellInfo$cellType)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| grnboost2 --num_workers | available cores - 1 | Parallel GRN inference |
| motif annotation source | Species-specific | Plant motifs differ from vertebrate |
| auc_threshold | 0.05 (per regulon) | Bin AUC binarization threshold |
| geneFiltering minCounts | 3 | Min UMI/count per gene |
| minGenes | 20 | Min genes per co-expression module |

## Plant-Specific Notes

- **Motif databases**: Standard cisTarget databases are vertebrate/human-specific.
  For plants, use PlantTFDB motifs converted to cisTarget format, or Arabidopsis
  motif collections (DAP-seq, ChIP-seq from Plant Cistrome Database).
- **Regulon interpretation in polyploids**: Homeologous TFs may share identical or
  diverged binding motifs. Regulons may merge (if motifs conserved) or split
  (if neo/subfunctionalization occurred). Compare regulon gene sets between homeologs.
- **Regulon specificity score (RSS)**: Essential for tissue/condition-specific analysis.
  In plants, cell-type and stress-specific regulons are the most informative.
- **Non-model species note**: If no species-specific motif database is available,
  use the Arabidopsis motif set as proxy. Validate key regulons via GO enrichment
  of target genes to confirm biological relevance.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No motifs annotated for X TFs" | Motif DB doesn't cover your species | Use Arabidopsis motifs; accept lower recovery rate |
| Empty regulons after ctx | Motifs too stringent or wrong DB | Reduce NES threshold (--nes_threshold 2.5) |
| Memory error in AUCell | Too many cells/samples | Downsample to 5000 cells; run in batches |
| "cisTarget db not found" | .feather files missing | Download from aertslab resources; check path |
| Low regulon count (<50) | Expression too flat or too few TFs | Check TF annotation completeness; verify data normalization |
