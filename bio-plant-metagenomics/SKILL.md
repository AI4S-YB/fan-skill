---
name: bio-plant-metagenomics
description: >
  Plant metagenomics — from raw reads to MAGs, taxonomic profiling,
  and functional annotation. Covers assembly (MEGAHIT/metaSPAdes),
  binning (MetaBAT2), quality control (CheckM), and plant microbiome analysis.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: MEGAHIT
workflow: true
---

# Plant Metagenomics Analysis

Metagenomic analysis tailored for plant-associated microbiomes — rhizosphere, phyllosphere, endosphere, and soil.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant metagenomics expert |
| **hybrid** | Matrix first → fall back to expert |

## Before Analysis

1. Check input: paired-end FASTQ files + sample metadata (compartment, host species, treatment)
2. Run `bash bio-plant-infra/scripts/check_env.sh` for tool availability
3. Estimate data volume: `seqkit stats *.fastq.gz`

## Analysis Flow

### Step 1: Quality Control and Preprocessing
Adapter trimming, quality filtering, host read removal. See `decision-matrix.yaml` > `preprocessing`.

### Step 2: Metagenome Assembly
Assemble reads into contigs. See `decision-matrix.yaml` > `assembly`.

### Step 3: Binning
Group contigs into metagenome-assembled genomes (MAGs). See `decision-matrix.yaml` > `binning`.

### Step 4: MAG Quality Assessment
Check MAG completeness and contamination. See `decision-matrix.yaml` > `mag_qc`.

### Step 5: Functional Annotation
Annotate MAGs and contigs with functional categories. See `decision-matrix.yaml` > `annotation`.

### Step 6: Visualization and Reporting
Community composition, functional profiles, PCoA. See `decision-matrix.yaml` > `visualization`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After QC | Reads retained | > 80% |
| After host removal | Host reads removed | < 5% host reads remaining |
| Assembly | N50 | > 1,000 bp (metagenome) |
| Assembly | Total assembly size | Within expected range for soil/rhizosphere |
| MAG quality (High) | Completeness > 90%, Contamination < 5% | High-quality draft |
| MAG quality (Medium) | Completeness > 50%, Contamination < 10% | Medium-quality draft |
| MAG count | MAGs per sample | 5-50 for soil, 1-10 for phyllosphere |

## References

- `bio-plant-infra/references/species-cheatsheet.md`
- `bio-plant-infra/references/plant-databases.md` — microbiome databases
- `bio-plant-infra/references/qc-thresholds.yaml` — metagenomics section
- `references/metagenomics-plant-special.md`
