---
name: bio-plant-rnaseq
description: >
  Plant RNA-seq differential expression analysis — from count matrix
  to DEGs, functional enrichment, and co-expression networks (WGCNA).
  Covers DESeq2, edgeR, limma-voom, Plant Reactome/GO/KEGG enrichment.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: DESeq2
workflow: true
---

# Plant RNA-seq Analysis

Transcriptome analysis tailored for plant species — differential expression, functional enrichment, co-expression.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant transcriptomics expert |
| **hybrid** | Matrix first → fall back to expert |

## Before Analysis

1. Check input: count matrix (genes × samples) + sample metadata (conditions/groups)
2. Run `bash bio-plant-infra/scripts/check_env.sh` for R package availability
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`

## Analysis Flow

### Step 1: Count Matrix Import
Detect format and import. See `decision-matrix.yaml` > `count_matrix_format`.

### Step 2: Low-Expression Filtering
Filter genes with insufficient counts. See `decision-matrix.yaml` > `low_expression_filter`.

### Step 3: Normalization + DE Analysis
Select tool via decision system. All methods in `tool-catalog/`.
Default: DESeq2 for ≥3 replicates.

### Step 4: Multiple Testing Correction
Benjamini-Hochberg FDR (α = 0.05).

### Step 5: Functional Enrichment
Plant Reactome → GO → KEGG, based on species availability.
See `tool-catalog/enrichment.md`.

### Step 6: Co-expression Network (Optional)
WGCNA for ≥15 samples with trait data.
See `tool-catalog/wgcna.md`.

### Step 7: Visualization + Report
Volcano plots, heatmaps, enrichment dot plots.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After filtering | Genes retained | > 10,000 |
| After DE | DEGs per contrast | ≥ 5 (may be 0 for weak effects) |
| After enrichment | Significant terms | ≥ 1 (may be 0 for few DEGs) |
| After WGCNA | Scale-free R² | > 0.80 |

## References

- `bio-plant-infra/references/species-cheatsheet.md`
- `bio-plant-infra/references/plant-databases.md` — expression DBs, enrichment tools
- `bio-plant-infra/references/qc-thresholds.yaml` — rnaseq section
- `references/rnaseq-plant-special.md`
