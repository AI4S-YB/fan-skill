---
name: bio-plant-hybrid-prediction
description: >
  Plant hybrid prediction — GCA/SCA estimation, heterotic group
  classification, genomic prediction of unobserved crosses, and
  optimal mate selection. Covers NCII and diallel designs, genomic
  BLUP, BayesR, and breeding constraints.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: ASReml-R / sommer
workflow: true
---

# Plant Hybrid Prediction

Hybrid breeding prediction for plant species — general and specific combining ability estimation, heterotic grouping, genomic cross prediction, and optimal mate selection with breeding constraints.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

| Mode | Behavior |
|------|----------|
| **rule** | Follow `decision-matrix.yaml` conditions → methods |
| **expert** | Read `analyst-notebook.md` and reason like a hybrid breeding expert |
| **hybrid** | Matrix first → fall back to expert |

## Before Analysis

1. Check input: phenotype data (yield, quality traits) + genotype data (SNP/DArT) + pedigree + experimental design
2. Run `bash bio-plant-infra/scripts/check_env.sh` for R package availability (sommer, BGLR, ASReml-R)
3. Verify mating design: NCII (factorial between heterotic groups), full/half diallel, or single-cross

## Analysis Flow

### Step 1: Preprocessing
Quality control of genotype and phenotype data. See `tool-catalog/preprocessing.md`.

### Step 2: GCA/SCA Estimation
Estimate general and specific combining ability from hybrid trials. See `decision-matrix.yaml` > `gca_sca`.

### Step 3: Heterotic Group Classification
Group inbred lines into heterotic pools using genomic or pedigree data. See `decision-matrix.yaml` > `heterotic_group`.

### Step 4: Cross Prediction
Predict performance of unobserved crosses using genomic prediction. See `decision-matrix.yaml` > `cross_prediction`.

### Step 5: Mate Selection
Optimally select parental combinations with breeding constraints. See `decision-matrix.yaml` > `mate_selection`.

### Step 6: Visualization and Reporting
Heterotic group PCA, prediction accuracy, genetic gain trajectories. See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| Genotype QC | Call rate | > 90% |
| Genotype QC | MAF | > 0.01 |
| Phenotype | Heritability (h2) | > 0.1 |
| GCA model | GCA/SCA variance ratio | GCA > SCA typical for yield |
| Cross prediction | Prediction accuracy (r) | > 0.3 (meaningful) |
| Heterotic groups | Silhouette score | > 0.3 |

## References

- `bio-plant-infra/references/species-cheatsheet.md`
- `bio-plant-infra/references/plant-databases.md` — breeding databases
- `bio-plant-infra/references/qc-thresholds.yaml` — hybrid prediction section
- `references/hybrid-prediction-plant-special.md`
