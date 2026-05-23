---
name: bio-plant-time-series
description: >
  Plant time-series expression analysis — from DESeq2-normalized counts
  through smoothing, clustering, and stage-specificity to publication-ready
  visualizations. Uses a dual-mode decision system (rule matrix + expert notebook).
tool_type: mixed
primary_tool: R
workflow: true
---

# Plant Time-Series Expression Analysis

Time-course gene expression analysis tailored for plant developmental
series, diurnal/circadian studies, and stress-response time courses.

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions -> methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant biologist |
| **hybrid** | Matrix first -> fall back to expert when no rule matches |

## Before Analysis

1. Profile the data: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`
4. Input should be DESeq2-normalized counts (variance-stabilizing transformation
   or regularized log) with time-point metadata

## Analysis Flow

### Step 1: Preprocessing
Smooth expression trajectories to separate signal from noise.
Refer to `decision-matrix.yaml` > `preprocessing` for method selection.
Methods in `tool-catalog/preprocessing.md`.

### Step 2: Clustering
Group genes with similar temporal expression patterns.
Refer to `decision-matrix.yaml` > `clustering` for method selection.
Methods in `tool-catalog/mfuzz.md` and `tool-catalog/kmeans-cluster.md`.

### Step 3: Trend Analysis
Identify statistically significant temporal trends (up, down, transient).
Refer to `decision-matrix.yaml` > `trend_analysis` for method selection.
Methods in `tool-catalog/trend-analysis.md`.

### Step 4: Stage Specificity
Quantify how specifically a gene is expressed in particular tissues or stages.
Refer to `decision-matrix.yaml` > `stage_specificity` for method selection.
Methods in `tool-catalog/stage-specificity.md`.

### Step 5: Visualization
Generate publication-ready figures.
Refer to `decision-matrix.yaml` > `visualization` for method selection.
Methods in `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After smoothing | Residual autocorrelation | < 0.3 |
| After clustering | Silhouette score | > 0.4 |
| After clustering | Cluster reproducibility (bootstrap) | > 0.8 |
| Tau index | Variance explained by top tissue | > 60% for tissue-specific genes |
| Circadian check | Periodogram peak significance | p < 0.05 if diurnal expected |

## References

- `bio-plant-infra/references/species-cheatsheet.md` -- species-specific guidance
- `bio-plant-infra/references/plant-databases.md` -- expression atlas sources
- `references/time-series-plant-special.md` -- plant developmental stages, tissue atlases, circadian handling
