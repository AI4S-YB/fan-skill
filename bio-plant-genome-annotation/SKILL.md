---
name: bio-plant-genome-annotation
description: >
  Plant genome annotation — from repeat masking to gene prediction
  and functional annotation. Covers RepeatMasker, BRAKER3, EggNOG-mapper,
  InterProScan, and BUSCO quality assessment.
  Uses a dual-mode decision system (rule matrix + expert notebook).
tool_type: mixed
primary_tool: repeatmasker
workflow: true
---

# Plant Genome Annotation

Comprehensive structural and functional annotation of plant genomes.

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant genomicist |
| **hybrid** | Matrix first → fall back to expert when no rule matches |

## Before Analysis

1. Profile the assembly: `bash bio-plant-infra/scripts/inspect_data.sh <genome_fasta>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`

## Analysis Flow

### Step 1: Preprocessing
Check genome assembly quality — contig N50, BUSCO completeness, GC content.
Refer to `tool-catalog/preprocessing.md` for assembly QC and masking preparation.

### Step 2: Repeat Masking
Mask repetitive elements using RepeatMasker with plant-specific repeat libraries.
Refer to `decision-matrix.yaml` > `repeat_masking`.

### Step 3: Gene Prediction
Predict protein-coding genes using BRAKER3 with evidence (RNA-seq or protein hints).
Refer to `decision-matrix.yaml` > `gene_prediction`.

### Step 4: Functional Annotation
Annotate predicted proteins with EggNOG-mapper and InterProScan.
Refer to `decision-matrix.yaml` > `functional_annotation`.

### Step 5: Quality Control
Assess annotation completeness with BUSCO in protein mode.
Refer to `decision-matrix.yaml` > `qc`.

### Step 6: Visualization
Generate annotation summary plots and genome browser tracks.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After repeat masking | Genome masked | 20-80% (depends on genome size) |
| After gene prediction | Gene count | Within expected range for species |
| After gene prediction | Average gene length | Comparable to related species |
| After functional annotation | Annotated proteins | > 70% with functional terms |
| After BUSCO | Complete BUSCOs | > 80% (embryophyta_odb10) |
| After BUSCO | Fragmented BUSCOs | < 10% |

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific guidance
- `bio-plant-infra/references/plant-databases.md` — genome & annotation sources
- `bio-plant-infra/references/qc-thresholds.yaml` — QC thresholds
- `references/plant-annotation-guide.md` — plant annotation specifics
