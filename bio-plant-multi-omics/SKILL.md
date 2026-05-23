---
name: bio-plant-multi-omics
description: >
  Plant multi-omics integration — transcriptomics, proteomics, metabolomics,
  epigenomics. Covers MOFA2 (unsupervised factor analysis), DIABLO (supervised
  classification), mixOmics PLS/sPLS (pairwise), data preprocessing, factor
  interpretation, and publication-quality visualization. Tailored for plant-specific
  challenges: polyploid subgenome handling, tissue/developmental stage integration,
  and non-model species strategies.
tool_type: mixed
primary_tool: mofa2
workflow: true
---

# Plant Multi-Omics Integration

Multi-omics integration for plant/crop species — from raw omics matrices to
latent factors, discriminative models, and biological interpretation.

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

1. Profile each omics data layer: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`
4. Review plant-specific multi-omics considerations: `references/multi-omics-plant-special.md`

## Analysis Flow

### Step 1: Data Preprocessing
- **Filter to common samples**: Only retain samples present in ALL omics layers (`tool-catalog/data-preprocessing.md`)
- **Normalize within each omics**: Scale/normalize each omics layer independently
- **Handle missing values**: Imputation strategies per omics type
- Refer to `decision-matrix.yaml` > `data_preprocessing`

### Step 2: Integration Method Selection
Select method via the decision system (`decision-matrix.yaml` > `integration_method`):

- **2 omics only**: mixOmics PLS/sPLS (`tool-catalog/mixomics.md`)
- **>= 3 omics, exploratory**: MOFA2 — learns latent factors (`tool-catalog/mofa.md`)
- **>= 2 omics, supervised classification**: DIABLO — maximizes class separation (`tool-catalog/diablo.md`)
- **No class labels**: MOFA2 or mixOmics PLS

### Step 3: Factor Interpretation
After MOFA2 factor discovery:
- Check variance explained by each factor (`tool-catalog/factor-analysis.md`)
- Identify top contributing features per factor
- Enrichment analysis for factor-associated genes/metabolites
- Refer to `decision-matrix.yaml` > `factor_interpretation`

### Step 4: Visualization + Report
- Multi-omics heatmaps and factor plots
- DIABLO sample plots and circos plots
- MOFA2 variance decomposition and factor-feature heatmaps
- See `tool-catalog/visualization.md`

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After preprocessing | Common samples retained | >= 50% of original |
| After preprocessing | Features per omics | >= 100 |
| After MOFA2 | Variance explained (top 5 factors) | >= 30% cumulatively |
| After DIABLO | Classification error rate | < 0.30 |
| After DIABLO | ncomp selected | >= 2 |
| After integration | Omics contributing to top factor | All omics layers present |

## Plant-Specific Workflow Notes

- **Polyploids**: Run subgenome-level MOFA2 separately. Check for homeolog coordination.
- **Non-model species**: Use homology-based annotation for feature interpretation.
- **Tissue/development**: Include tissue or developmental stage as covariates if available.
- **Multi-environment**: If samples span environments, treat environment as a DIABLO outcome or MOFA2 covariate.

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific guidance
- `bio-plant-infra/references/plant-databases.md` — genome & annotation sources
- `bio-plant-infra/references/qc-thresholds.yaml` — QC thresholds
- `references/multi-omics-plant-special.md` — plant multi-omics specifics
