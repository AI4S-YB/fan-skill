---
name: bio-plant-genomic-selection
description: >
  Plant genomic selection / genomic prediction — from genotype+phenotype
  to GEBVs and selection decisions. Covers rrBLUP, GBLUP, BayesA/B/C,
  deep learning GS, cross-validation, and multi-trait indices.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: sommer
workflow: true
---

# Plant Genomic Selection

Genomic prediction and selection for plant breeding — estimate breeding values, cross-validate, rank selections.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Follow `decision-matrix.yaml` conditions → models |
| **expert** | Read `analyst-notebook.md` and reason like a plant breeder |
| **hybrid** | Matrix first → fall back to expert |

## Before Analysis

1. Profile genotype data: `bash bio-plant-infra/scripts/inspect_data.sh <input>`
2. Check phenotype data format (sample ID + trait values)
3. Estimate trait heritability (if possible from replicates or literature)
4. Run `bash bio-plant-infra/scripts/check_env.sh`

## Analysis Flow

### Step 1: Genotype Preprocessing
QC filtering, MAF, missing rate. See `decision-matrix.yaml` > `genotype_preprocessing`.

### Step 2: Model Selection
Choose GS model via decision system.
See `decision-matrix.yaml` > `model_selection`.

### Step 3: Cross-Validation
Estimate prediction accuracy. See `decision-matrix.yaml` > `cross_validation`.
**Critical**: Use stratified CV for structured populations.

### Step 4: Accuracy Assessment
Pearson r + RMSE. See `decision-matrix.yaml` > `accuracy_metrics`.

### Step 5: Breeding Value Ranking
Rank genotypes by GEBV, apply selection strategy.
See `decision-matrix.yaml` > `selection_strategy`.

### Step 6: Visualization + Report
Accuracy boxplots, predicted vs observed, GEBV distribution.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After preprocessing | Markers retained | > 500 |
| After CV | Pearson r | > 0.30 |
| After CV | CV fold SD | < 0.15 |

## References

- `bio-plant-infra/references/species-cheatsheet.md`
- `bio-plant-infra/references/qc-thresholds.yaml` — genomic_selection section
- `references/gs-plant-special.md`
