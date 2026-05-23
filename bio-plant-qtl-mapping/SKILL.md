---
name: bio-plant-qtl-mapping
description: >
  Plant QTL mapping for biparental populations — genetic map construction,
  composite interval mapping (CIM), multiple QTL mapping (MQM),
  and multi-environment QTL analysis. Dual-mode B+C decision system.
tool_type: mixed
primary_tool: R/qtl
workflow: true
---

# Plant QTL Mapping

Linkage-based QTL mapping for biparental populations (F2, RIL, DH, BC).

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant geneticist |
| **hybrid** | Matrix first → fall back to expert when no rule matches |

## Before Analysis

1. Profile the data: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. Verify population type (F2/RIL/DH/BC) and marker encoding
4. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`

## Analysis Flow

### Step 1: Genetic Map Construction
- Import genotype data into R/qtl (`read.cross`)
- Check marker segregation distortion
- Estimate genetic map distance (`est.map`) or use ASMap / LepMap3 for large datasets
- Refer to `decision-matrix.yaml` > `genetic_map`

### Step 2: QTL Scanning
- Select method: IM / CIM / MQM
- Refer to `decision-matrix.yaml` > `qtl_method`
- All methods in `tool-catalog/` with full parameters and common errors

### Step 3: LOD Threshold Determination
- Permutation test (1000 permutations for n >= 100)
- Analytical threshold for small samples
- Refer to `decision-matrix.yaml` > `lod_threshold`

### Step 4: Multi-Environment QTL
- Joint analysis or per-environment scanning
- Refer to `decision-matrix.yaml` > `multi_env_qtl`

### Step 5: QTL Interval Annotation and Visualization
- LOD profile plots, genetic map plots, QTL effect plots
- See `tool-catalog/visualization.md`

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| Genetic map | Marker retention | > 80% |
| Genetic map | Segregation distortion | < 20% markers |
| Genetic map | Map length per chromosome | < 300 cM |
| QTL scan | LOD threshold | Permutation-based |
| QTL scan | QTL confidence interval | < 30 cM |

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific guidance
- `bio-plant-infra/references/plant-databases.md` — genome & annotation sources
- `references/qtl-plant-special.md` — plant QTL specifics
