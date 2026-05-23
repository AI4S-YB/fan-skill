---
name: bio-plant-metabolomics
description: >
  Plant metabolomics analysis — from LC-MS/GC-MS peak detection
  to differential metabolite analysis and pathway mapping.
  Covers XCMS, MZmine, limma, SIRIUS+CSI:FingerID, and
  PlantCyc/KEGG pathway enrichment.
  Uses a dual-mode decision system (rule matrix + expert notebook).
tool_type: mixed
primary_tool: xcms
workflow: true
---

# Plant Metabolomics Analysis

Untargeted and targeted metabolomics for plant systems — secondary metabolites, primary metabolism, and specialized pathways.

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant metabolomist |
| **hybrid** | Matrix first → fall back to expert when no rule matches |

## Before Analysis

1. Profile the data: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. Identify platform type: LC-MS, GC-MS, or NMR

## Analysis Flow

### Step 1: Preprocessing
Raw file conversion (vendor format → mzML), quality check, and sample grouping.
Refer to `tool-catalog/preprocessing.md`.

### Step 2: Peak Detection and Alignment
Detect chromatographic peaks and align features across samples.
Refer to `decision-matrix.yaml` > `peak_detection`.

### Step 3: Differential Analysis
Identify metabolites significantly different between conditions.
Refer to `decision-matrix.yaml` > `diff_analysis`.

### Step 4: Metabolite Annotation
Putative and confident identification of metabolite features.
Refer to `decision-matrix.yaml` > `annotation`.

### Step 5: Pathway Mapping
Map significant metabolites to plant metabolic pathways.
Refer to `decision-matrix.yaml` > `pathway_mapping`.

### Step 6: Visualization
Generate PCA scores plots, volcano plots, heatmaps, and pathway diagrams.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After preprocessing | TIC variation (QC samples) | CV < 30% |
| After peak detection | Features detected | > 1000 (LC-MS untargeted) |
| After alignment | Retention time drift | < 5s (within batch) |
| After differential analysis | Features with p < 0.05 | > 0 (can be low for subtle phenotypes) |
| After pathway mapping | Features mapped to pathways | > 10% of significant features |

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific metabolism notes
- `bio-plant-infra/references/plant-databases.md` — metabolomics databases
- `bio-plant-infra/references/qc-thresholds.yaml` — QC thresholds
- `references/plant-metabolomics-guide.md` — plant metabolomics specifics
