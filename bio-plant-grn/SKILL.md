---
name: bio-plant-grn
description: >
  Plant gene regulatory network inference — from expression matrix
  to TF annotation, network reconstruction (GENIE3/SCENIC), hub gene
  detection, module analysis, and visualization. Tailored for plant-specific
  TF families, polyploid homeolog handling, and non-model species strategies.
tool_type: mixed
primary_tool: genie3
workflow: true
---

# Plant Gene Regulatory Network (GRN) Inference

Gene regulatory network inference for plant/crop species — expression matrix to
functional modules.

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions -> methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant systems biologist |
| **hybrid** | Matrix first -> fall back to expert when no rule matches |

## Before Analysis

1. Profile the expression data: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`
4. Review plant TF specifics: `references/grn-plant-special.md`

## Analysis Flow

### Step 1: Expression Data Preparation
- Input: normalized count matrix (genes x samples) in TSV/CSV
- Filter lowly expressed genes (mean count >= 1, or TPM >= 1 in >= 20% of samples)
- Log-transform if needed (log2(TPM+1) or variance-stabilizing transformation)
- Check for batch effects — PCA on expression matrix before and after correction

### Step 2: Transcription Factor Annotation
- Identify which genes are TFs in your species
- Refer to `decision-matrix.yaml` > `tf_database` for database selection
- For model species: PlantTFDB direct lookup
- For non-model species: homology-based inference from Arabidopsis
- Supporting tool: `tool-catalog/tf-database.md`

### Step 3: Network Inference Method Selection
Select method via the decision system (`decision-matrix.yaml` > `inference_method`):
- **GENIE3**: Tree-based, >= 15 samples, TF annotation available. Directed edges TF->target.
- **SCENIC**: >= 20 samples, adds regulon detection and activity scoring.
- **WGCNA**: Undirected co-expression, no TF annotation needed.
All methods documented in `tool-catalog/`.

### Step 4: Network Analysis
After network construction:
- **Hub gene detection**: degree + betweenness centrality (`tool-catalog/hub-gene.md`)
- **Module detection**: community detection for networks >= 500 nodes (`tool-catalog/network-analysis.md`)
- Refer to `decision-matrix.yaml` > `network_analysis`

### Step 5: Visualization
- TF-target subnetwork plots
- Hub gene rankings
- Module structure visualization
- See `tool-catalog/visualization.md`

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After filtering | Genes retained | >= 5000 |
| After filtering | Samples retained | 100% |
| TF annotation | TFs identified | >= 100 |
| Network edge count | Edges | >= 1000 |
| Hub detection | Top 50 hub TFs | >= 1 known TF family hit |
| Module detection | Modules found | >= 2 |

## Plant-Specific Workflow Notes

- **Polyploids**: Run TF annotation per subgenome. Check for homeolog pairs in hub gene lists.
- **Non-model species**: Use homology-based TF inference. Validate top edges against known interactions.
- **TF family enrichment**: Plant-specific families (MYB, WRKY, NAC, bHLH, AP2/ERF, GRAS, MADS-box) should dominate top hubs.

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific guidance
- `bio-plant-infra/references/plant-databases.md` — genome & annotation sources
- `bio-plant-infra/references/qc-thresholds.yaml` — QC thresholds
- `references/grn-plant-special.md` — plant GRN specifics (TF families, polyploidy)
