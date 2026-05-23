---
name: bio-plant-pan-genome
description: >
  Plant pan-genome analysis — graph construction (Minigraph/PGGB),
  core/variable gene classification, PAV detection, SV genotyping,
  and pan-genome visualization. Supports small to large plant genome cohorts.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: minigraph
workflow: true
---

# Plant Pan-Genome Analysis

Pan-genome analysis for plant species — graph construction, presence/absence variation, structural variant genotyping, and comparative visualization.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant pan-genomics expert |
| **hybrid** | Matrix first → fall back to expert |

## Before Analysis

1. Check input: genome FASTA list + annotation GFF list + sample metadata
2. Run `bash bio-plant-infra/scripts/check_env.sh` for tool availability
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`

## Analysis Flow

### Step 1: Pan-genome Graph Construction
Build a sequence graph from multiple plant genome assemblies. See `decision-matrix.yaml` > `graph_construction`.

### Step 2: Core/Variable Gene Classification
Classify genes by presence frequency across genomes. See `decision-matrix.yaml` > `core_variable`.

### Step 3: PAV Detection
Detect presence/absence variation across the pan-genome. See `decision-matrix.yaml` > `pav_detection`.

### Step 4: SV Genotyping
Genotype structural variants in new samples by mapping to the pan-genome graph. See `decision-matrix.yaml` > `sv_genotyping`.

### Step 5: Visualization
Render pan-genome graphs and PAV distribution. See `decision-matrix.yaml` > `visualization`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| Graph construction | Nodes in graph | > 1,000,000 expected for plant genomes |
| Core genome | Core gene count | Typically 30-60% of total gene families |
| PAV matrix | Genes with PAV | Usually 15-40% of pangenome shows PAV |
| SV genotyping | Genotyping rate | > 85% |
| Graph visualization | Graph connectivity | Reasonable graph structure, no disconnected components |

## References

- `bio-plant-infra/references/species-cheatsheet.md`
- `bio-plant-infra/references/plant-databases.md` — genome databases
- `bio-plant-infra/references/qc-thresholds.yaml` — pan-genome section
- `references/pan-genome-plant-special.md`
