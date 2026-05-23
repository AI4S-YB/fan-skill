---
name: bio-plant-proteomics
description: >
  Plant proteomics analysis — from raw MS spectra to differential
  protein abundance, PPI networks, and PTM analysis. Covers DDA
  (MaxQuant) and DIA (DIA-NN) workflows, limma-based differential
  analysis, STRING plant PPI, and phosphorylation analysis.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: MaxQuant
workflow: true
---

# Plant Proteomics Analysis

Proteomics analysis tailored for plant species — label-free quantification, differential abundance, protein-protein interaction networks, and post-translational modification analysis.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant proteomics expert |
| **hybrid** | Matrix first → fall back to expert |

## Before Analysis

1. Check input: raw MS files (.raw, .mzML, .d) + sample metadata + experimental design
2. Verify acquisition mode: DDA (data-dependent acquisition) or DIA (data-independent acquisition)
3. Check protein database availability for your plant species

## Analysis Flow

### Step 1: Protein Identification and Quantification
Search raw spectra against plant protein database. See `decision-matrix.yaml` > `quantification`.

### Step 2: Data Preprocessing
Normalization, missing value imputation, log2 transformation. See `tool-catalog/diff-analysis.md`.

### Step 3: Differential Abundance Analysis
Identify significantly changed proteins between conditions. See `decision-matrix.yaml` > `diff_analysis`.

### Step 4: Protein-Protein Interaction Network
Map differentially abundant proteins to plant PPI networks. See `decision-matrix.yaml` > `ppi_network`.

### Step 5: PTM Analysis (Optional)
Analyze phosphorylation or other PTMs. See `decision-matrix.yaml` > `ptm`.

### Step 6: Visualization and Reporting
Volcano plots, heatmaps, PPI networks. See `decision-matrix.yaml` > `visualization`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| Peptide identification | Peptides identified | > 10,000 (plant proteome) |
| Protein identification | Protein groups | > 2,000 |
| Missing values | Proteins with >50% valid values | < 20% per condition |
| PCA | Clustering by condition | Expected grouping visible |
| DEPs | Significant proteins (FDR < 0.05) | >= 5 (may be 0 for weak effects) |

## References

- `bio-plant-infra/references/species-cheatsheet.md`
- `bio-plant-infra/references/plant-databases.md` — protein databases
- `bio-plant-infra/references/qc-thresholds.yaml` — proteomics section
- `references/proteomics-plant-special.md`
