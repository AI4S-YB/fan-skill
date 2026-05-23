---
name: bio-plant-eqtl
description: >
  Plant expression quantitative trait locus (eQTL) analysis — from genotype
  + expression data to cis/trans-eQTL identification, colocalization with GWAS
  loci, and multi-tissue eQTL comparison. Uses MatrixEQTL for fast linear model
  testing, with plant-specific guidance on polyploid mapping complexity and
  tissue atlas availability.
tool_type: mixed
primary_tool: matrixeqtl
workflow: true
---

# Plant eQTL Analysis

Expression quantitative trait locus analysis tailored for plant/crop species —
linking genetic variants to gene expression variation.

## Decision Modes

Set in `params.yaml`:

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Strictly follow `decision-matrix.yaml` conditions -> methods |
| **expert** | Read `analyst-notebook.md` and reason like a plant geneticist |
| **hybrid** | Matrix first -> fall back to expert when no rule matches |

## Before Analysis

1. Profile the data: `bash bio-plant-infra/scripts/inspect_data.sh <input_path>`
2. Check the environment: `bash bio-plant-infra/scripts/check_env.sh`
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`
4. Review plant eQTL specifics: `references/eqtl-plant-special.md`

## Analysis Flow

### Step 1: Data Preparation and QC

- Genotype data: VCF or PLINK binary format. Filter MAF >= 0.05, call rate >= 0.9.
- Expression data: normalized count matrix (genes x samples) in TSV/CSV.
  Filter lowly expressed genes (>= 1 TPM in >= 20% of samples).
- Match sample IDs between genotype and expression datasets.
- Account for hidden covariates (PEER factors, PCA on expression).
- See `tool-catalog/preprocessing.md` for detailed R workflow.

### Step 2: eQTL Mapping Method Selection

Select method via the decision system (`decision-matrix.yaml` > `eqtl_type`):
- **cis-eQTL**: Variants within +/- 1 Mb of gene TSS. Use MatrixEQTL linear model.
- **trans-eQTL**: Variants on different chromosomes or > 1 Mb from TSS.
  Requires large sample size (>= 200).
- All methods documented in `tool-catalog/`.

### Step 3: Multiple Testing Correction

Millions of tests demand rigorous correction:
- **FDR < 0.05** as primary threshold for cis-eQTL.
- **Permutation-based FDR** recommended for >= 1M tests to control false positives.
- Refer to `decision-matrix.yaml` > `multiple_testing`.
- See `tool-catalog/cis-trans.md` for full threshold guidance.

### Step 4: Colocalization with GWAS

If GWAS summary statistics are available:
- **coloc**: Bayesian colocalization testing (H0-H4 posterior probabilities).
- Useful for linking eQTL signals to known GWAS loci.
- See `decision-matrix.yaml` > `colocalization` and `tool-catalog/colocalization.md`.

### Step 5: Multi-Tissue eQTL Analysis

When expression data from >= 2 tissues is available:
- Classify eQTLs as **tissue-specific** or **tissue-shared**.
- Use multivariate adaptive shrinkage (mashR) or meta-analysis approaches.
- See `decision-matrix.yaml` > `multi_tissue` and `tool-catalog/multi-tissue.md`.

### Step 6: Visualization

- eQTL Manhattan plot per gene / per chromosome.
- cis vs trans classification summary.
- Colocalization posterior probability bar plots.
- Tissue-sharing heatmaps and UpSet plots.
- See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After genotype QC | Variant retention | >= 70% |
| After expression QC | Genes retained | >= 5000 |
| After expression QC | Samples with both data types | >= 80% |
| After eQTL mapping | cis-eGenes (FDR < 0.05) | >= 0 (0 can be valid for small n) |
| After eQTL mapping | Median cis-eQTL distance | < 50 kb |
| Colocalization | PP.H4 >= 0.75 for shared signals | >= 1 locus (if GWAS overlap) |

## Plant-Specific Workflow Notes

- **Polyploids**: Map eQTL per subgenome. Homeolog expression bias can confound trans-eQTL detection.
  Use subgenome-specific SNP annotation to avoid cross-homeolog mapping.
- **Tissue atlases**: Leverage species-specific expression atlases (e.g., Rice eFP, Maize Atlas,
  SoyBase, Wheat EXPVIZ) for tissue-specific eQTL interpretation.
- **Low sample size**: Plant eQTL studies often have limited samples (n < 100). cis-eQTL
  detection power is reasonable; trans-eQTL requires much larger cohorts.
- **Non-model species**: When no reference eQTL dataset exists, compare to Arabidopsis/rice
  orthologs for functional validation.

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific guidance
- `bio-plant-infra/references/plant-databases.md` — genome & annotation sources
- `bio-plant-infra/references/qc-thresholds.yaml` — QC thresholds
- `references/eqtl-plant-special.md` — plant eQTL specifics (polyploidy, tissue atlases)
