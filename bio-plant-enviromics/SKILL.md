---
name: bio-plant-enviromics
description: >
  Plant enviromics analysis — environmental characterization, G×E modeling,
  TPE (target population of environments) delineation, and genotype stability.
  Covers PCA-based environmental indices, AMMI/GGE, factor-analytic models,
  reaction norms, and multi-metric stability analysis.
  Uses a dual-mode decision system (rule matrix + expert notebook).
tool_type: mixed
primary_tool: R
workflow: true
---

# Plant Enviromics Analysis

Environmental characterization and genotype-by-environment interaction (G×E)
analysis tailored for plant breeding and crop science.

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant enviromics scientist |
| **hybrid** | Matrix first → fall back to expert when no rule matches |

## Before Analysis

1. Profile the data: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`

## Analysis Flow

### Step 1: Environmental Index Construction
Characterize environments using climate/soil variables.
Refer to `decision-matrix.yaml` > `env_index` for method selection.
When ≥ 3 climate variables exist → PCA-based composite index.
When < 3 variables → use raw variables directly.
See `tool-catalog/env-index.md`.

### Step 2: Climate Data Integration
Integrate external climate/weather data with trial metadata.
Handles gridded climate products, station data, and soil grids.
See `tool-catalog/climate-integration.md`.

### Step 3: G×E Modeling
Select method via decision system based on environment count and genotype count.
- ≥ 5 env + ≥ 50 genotypes → Factor-Analytic (FA) models
- ≥ 2 env → AMMI/GGE biplot analysis
- Environmental index available → Reaction norm models
All methods in `tool-catalog/gxe-modeling.md`.

### Step 4: TPE Analysis
When ≥ 5 environments, cluster environments into Target Population of
Environments (TPE) groups for breeding strategy design.
See `tool-catalog/tpe-analysis.md`.

### Step 5: Stability Analysis
Compute multiple stability metrics for cross-validation.
See `tool-catalog/stability.md`.

### Step 6: Visualization + Report
Biplots, heatmaps, stability plots, environmental maps.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After env-index | PCA variance explained (PC1+PC2) | > 60% |
| After env-index | Climate variable coverage | > 80% |
| After G×E | GE variance ratio | > 5% of total |
| After G×E | Model convergence | TRUE |
| After stability | Metric agreement (Kendall W) | > 0.6 |
| After TPE | Silhouette score | > 0.3 |

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific guidance
- `bio-plant-infra/references/plant-databases.md` — climate & soil data sources
- `references/enviromics-plant-special.md` — plant enviromics specifics
