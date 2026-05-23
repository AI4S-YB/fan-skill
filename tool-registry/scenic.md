# SCENIC — Single-Cell Regulatory Network Inference and Clustering

**Goal:** From expression matrix to regulons — infers co-expression modules between
TFs and targets (via GENIE3/GRNBoost2), prunes indirect targets via cis-regulatory
motif enrichment (RcisTarget), and scores regulon activity per sample (AUCell).

**Best for:** >= 20 samples, TF annotation available, regulon activity needed.

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
