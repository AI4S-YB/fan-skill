---
name: bio-plant-marker
description: >
  Plant molecular marker development — from GWAS peaks/QTL to PCR markers.
  Covers KASP, InDel, SSR, and functional (dCAPS/CAPS) marker design,
  validation, and parental combination recommendation.
  Dual-mode B+C decision system.
tool_type: mixed
primary_tool: Primer3
workflow: true
---

# Plant Molecular Marker Development

From variant to validated breeding marker — design, validate, deploy.

## Decision Modes

```yaml
decision_mode: hybrid
```

## Analysis Flow

1. **Variant Type Detection**: SNP, InDel, SSR, or functional variant?
2. **Marker Type Selection**: KASP, gel-based, or functional marker?
3. **Primer Design**: Primer3 with plant-specific parameters
4. **Marker Validation**: In-silico PCR + polymorphism check
5. **Parental Recommendation**: Complementary alleles, GEBV pairing, diversity maximization

## QC Checkpoints

| Checkpoint | Metric | Pass |
|------------|--------|------|
| After primer design | Primer pairs found | >= 1 per target |
| After validation | Specific amplification | Yes (single locus) |
| After polymorphism | MAF in population | > 0.05 |
