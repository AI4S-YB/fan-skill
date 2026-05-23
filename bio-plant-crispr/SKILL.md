---
name: bio-plant-crispr
description: >
  Plant CRISPR/Cas9 design — sgRNA design (CRISPOR/CRISPR-P),
  off-target prediction (Cas-OFFinder), editing efficiency
  prediction (DeepSpCas9), HDR template design, and validation
  strategy. Covers major crops and custom genomes.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: CRISPOR
workflow: true
---

# Plant CRISPR/Cas9 Design

CRISPR/Cas9 genome editing design for plant species -- sgRNA design, off-target analysis, editing efficiency prediction, HDR template design, and validation.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant genome editing expert |
| **hybrid** | Matrix first → fall back to expert |

## Before Analysis

1. Check input: target gene sequence (FASTA) + species + editing goal (KO, knock-in, base editing, prime editing)
2. Run `bash bio-plant-infra/scripts/check_env.sh` for tool availability
3. Verify Cas variant: SpCas9 (NGG PAM), Cas12a (TTTV PAM), or other

## Analysis Flow

### Step 1: sgRNA Design
Design guide RNAs targeting your gene of interest. See `decision-matrix.yaml` > `sgrna_design`.

### Step 2: Off-Target Analysis
Screen candidate sgRNAs against the whole genome for potential off-target sites. See `decision-matrix.yaml` > `off_target`.

### Step 3: Editing Efficiency Prediction
Predict on-target editing efficiency of each sgRNA. See `decision-matrix.yaml` > `editing_efficiency`.

### Step 4: HDR Template Design (Optional)
Design donor template for homology-directed repair knock-in. See `decision-matrix.yaml` > `hdr_template`.

### Step 5: Validation Strategy
Plan genotyping assays to confirm edits. See `decision-matrix.yaml` > `validation`.

### Step 6: Visualization
Generate sgRNA target maps and vector diagrams. See `task-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| sgRNA design | sgRNAs per gene | 2-5 candidates |
| Off-target | Off-target sites (0-3 mismatches) | < 10 total (preferably 0-1) |
| Efficiency | Predicted efficiency | > 50% (DeepSpCas9 score) |
| sgRNA specificity | CFD score | > 80 |
| PAM | Canonical PAM present | NGG (SpCas9) or TTTV (Cas12a) |
| GC content | sgRNA GC% | 40-60% |

## References

- `bio-plant-infra/references/species-cheatsheet.md`
- `bio-plant-infra/references/plant-databases.md` — genome databases for off-target search
- `bio-plant-infra/references/qc-thresholds.yaml` — CRISPR section
- `references/crispr-plant-special.md`
