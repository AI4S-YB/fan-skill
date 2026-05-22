---
name: bio-plant-comparative
description: >
  Plant comparative genomics — synteny (MCScanX), Ks analysis,
  gene family identification (OrthoFinder/HMM), selection pressure
  (PAML dN/dS), and species/gene tree reconciliation.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: MCScanX
workflow: true
---

# Plant Comparative Genomics

Comparative genomics for plant species — synteny, gene families, selection, evolution.

## Decision Modes

```yaml
decision_mode: hybrid
```

## Analysis Flow

1. **Synteny Analysis**: MCScanX intra/inter-species collinearity
2. **Ks Analysis**: Dating duplications via synonymous substitutions
3. **Gene Family**: OrthoFinder or HMM-based identification
4. **Selection Pressure**: PAML branch-site/site models for dN/dS
5. **Tree Reconciliation**: Species tree vs gene tree

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After synteny | Collinear blocks detected | ≥ 5 |
| After Ks | Peak identified | ≥ 1 |
| After gene family | Families identified | ≥ 100 |
| After selection | Genes with dN/dS > 1 | ≥ 0 |
