---
name: bio-plant-phenotype
description: >
  Plant phenotype data analysis — heritability estimation, BLUP/BLUE,
  spatial correction, multi-environment trial summary (AMMI/GGE).
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: sommer
workflow_status: active
---

# Plant Phenotype Analysis

Phenotype data analysis for plant breeding — from raw trial data to heritability estimates and genotype rankings.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid (recommended)
```

## Analysis Flow

1. **Descriptive Statistics**: Distribution, CV, outliers
2. **Heritability**: Mixed model (multi-env) or ANOVA (single-env)
3. **BLUP/BLUE**: Genotype effect estimation
4. **Spatial Correction**: AR1×AR1 or spline-based if field trial
5. **MET Summary**: AMMI/GGE for multi-environment trials

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| Data quality | Trait CV | < 30% |
| Heritability | H⊃2; or h⊃2; | > 0.10 for GS, > 0.20 for GWAS |
| Spatial | CV after correction | Lower than before |
| MET | Environments with data | ≥ 2 for cross-env analysis |
