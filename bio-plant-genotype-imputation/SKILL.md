---
name: bio-plant-genotype-imputation
description: >
  Plant genotype imputation — from post-QC genotype data
  to fully imputed variant set. Covers pre-imputation QC,
  Beagle 5, IMPUTE2, Minimac4, polyploid imputation, and
  post-imputation filtering. Uses a dual-mode decision system
  (rule matrix + expert notebook).
tool_type: mixed
primary_tool: beagle
workflow: true
---

# Plant Genotype Imputation

Genotype imputation tailored for plant/crop species.

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

## Analysis Flow

### Step 1: Pre-Imputation QC
Standard PLINK QC: MAF, call rate, HWE.
Refer to `decision-matrix.yaml` > `pre_imputation_qc` for parameter selection.
Check species type -- skip HWE for inbred species.
See `tool-catalog/pre-imputation-qc.md` for full commands.

### Step 2: Germplasm & Population Structure Assessment
- Check population structure: PCA via PLINK
- For structured populations: stratify before imputation
- For outcross species: ensure diverse reference panel
Refer to `decision-matrix.yaml` > `germplasm_consideration`.

### Step 3: Tool Selection
Select imputation tool via decision system.
All methods in `tool-catalog/` with full parameters and common errors.
Key decision factors: sample count, SNP count, reference panel availability, ploidy.

### Step 4: Phasing & Imputation
- Beagle: integrated phasing + imputation
- IMPUTE2: requires SHAPEIT pre-phasing + reference haplotype panel
- Minimac4: fastest for large samples, requires M3VCF reference
See `decision-matrix.yaml` > `phasing_strategy`.

### Step 5: Post-Imputation QC
Filter by R-squared (info score) and MAF.
Refer to `decision-matrix.yaml` > `post_imputation_qc` for parameter selection.
Tighten thresholds for GWAS downstream analysis.
See `tool-catalog/post-imputation-qc.md` for full commands.

### Step 6: Visualization & Report
Imputation accuracy plots, before/after density comparison, missing rate heatmap.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After pre-QC | Sample retention | > 90% |
| After pre-QC | Variant retention | > 50% |
| After imputation | Mean R² (info score) | > 0.7 |
| After imputation | Variants with R² > 0.3 | > 80% |
| After post-QC | Variants retained | > 60% of imputed |
| After post-QC | Concordance (if validation) | > 95% |

## References

- `bio-plant-infra/references/species-cheatsheet.md` -- species-specific guidance
- `bio-plant-infra/references/plant-databases.md` -- genome & annotation sources
- `bio-plant-infra/references/qc-thresholds.yaml` -- QC thresholds
- `references/imputation-plant-special.md` -- plant imputation specifics
