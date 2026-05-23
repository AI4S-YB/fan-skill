---
name: bio-plant-amplicon
description: >
  Plant amplicon sequencing analysis — 16S rRNA, ITS, and other marker
  gene surveys for plant-associated microbiomes. Covers DADA2/QIIME2
  denoising, taxonomy assignment (SILVA/UNITE), diversity analysis,
  and differential abundance testing (ANCOM-BC).
  Uses a dual-mode decision system (rule matrix + expert notebook).
tool_type: mixed
primary_tool: qiime2
workflow: true
---

# Plant Amplicon Sequencing Analysis

Amplicon-based microbiome profiling for plant-associated bacterial and fungal communities.

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a microbial ecologist |
| **hybrid** | Matrix first → fall back to expert when no rule matches |

## Before Analysis

1. Profile the data: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. If metadata available, read experimental design from `metadata.tsv`

## Analysis Flow

### Step 1: Denoising and ASV/OTU Table
Generate amplicon sequence variants (ASVs) or OTUs from raw reads.
Refer to `decision-matrix.yaml` > `denoising`.

### Step 2: Taxonomy Assignment
Assign taxonomy to representative sequences.
Refer to `decision-matrix.yaml` > `taxonomy`.

### Step 3: Diversity Analysis
Calculate alpha and beta diversity metrics.
Refer to `decision-matrix.yaml` > `diversity`.

### Step 4: Differential Abundance
Identify taxa differentially abundant between groups.
Refer to `decision-matrix.yaml` > `diff_abundance`.

### Step 5: Visualization
Generate publication-quality plots: bar plots, PCoA, heatmaps.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After denoising | Reads retained | > 50% |
| After denoising | ASVs detected | > 100 |
| After taxonomy | Classified to phylum | > 80% |
| After taxonomy | Classified to genus | > 30% (16S) or > 20% (ITS) |
| After diversity | Sequencing depth per sample | > 5000 reads |
| After diversity | Sample retention after rarefaction | > 90% |

## References

- `bio-plant-infra/references/species-cheatsheet.md` — plant host species guidance
- `bio-plant-infra/references/plant-databases.md` — plant microbiome databases
- `bio-plant-infra/references/qc-thresholds.yaml` — QC thresholds
- `references/plant-microbiome-guide.md` — plant amplicon specifics
