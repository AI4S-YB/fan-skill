---
name: bio-plant-gwas
description: >
  Plant genome-wide association study — from genotype + phenotype
  to candidate genes. Covers PLINK QC, population structure correction,
  GAPIT (CMLM/FarmCPU), BLINK, and MET-GWAS for multi-environment trials.
  Uses a dual-mode decision system (rule matrix + expert notebook).
tool_type: mixed
primary_tool: plink2
workflow: true
---

# Plant GWAS Analysis

Genome-wide association study tailored for plant/crop species.

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
3. If species detected, read `bio-plant-infra/references/species-cheatsheet.md`

## Analysis Flow

### Step 1: Data Conversion
- BED → skip
- VCF → `plink2 --vcf in.vcf.gz --make-bed --out prefix`
- Hapmap → GAPIT reads directly

### Step 2: Genotype QC
Standard PLINK QC: MAF, call rate, HWE.
Refer to `decision-matrix.yaml` > `qc_strategy` for parameter selection.
Check species type — skip HWE for inbred species.

### Step 3: Population Structure
Run PCA via PLINK. Check λ.
Refer to `decision-matrix.yaml` > `covariate_strategy`.

### Step 4: GWAS
Select method via decision system.
All methods in `tool-catalog/` with full parameters and common errors.

### Step 5: Visualization + Report
Manhattan, QQ plot, significant SNP table.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After QC | Sample retention | > 90% |
| After QC | Variant retention | > 50% |
| After PCA | PC1 variance | > 5% |
| After GWAS | λ | < 1.2 (inbred: < 1.3) |
| After GWAS | Significant SNPs | ≥ 0 (0 can be valid) |

## References

- `bio-plant-infra/references/species-cheatsheet.md` — species-specific guidance
- `bio-plant-infra/references/plant-databases.md` — genome & annotation sources
- `bio-plant-infra/references/qc-thresholds.yaml` — QC thresholds
- `references/gwas-plant-special.md` — plant GWAS specifics
