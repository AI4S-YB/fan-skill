# Ks Analysis (Synonymous Substitution Rate)

**Goal:** Estimate divergence time between gene pairs via synonymous substitutions
**Best for:** Dating whole-genome duplications, estimating gene pair ages

## Prerequisites
- KaKs_Calculator 2.0+ or PAML (codeml)
- Codon-aligned CDS sequences for each gene pair
- ParaAT.pl for batch processing

## Basic Usage

```bash
KaKs_Calculator -i aligned.cds.axt -o output.kaks -m YN
```

## Plant-Specific Notes
- Most plant genomes have 1-3 ancient WGD events
- Ks peak around 1.0-2.0 = ancient WGD shared by all eudicots
- Ks peak at 0.1-0.5 = recent lineage-specific WGD
- Grass-specific ρ WGD: Ks ≈ 0.8-1.2

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "Ks = NA" | Saturation (Ks > 2) | Flag as "saturated — age unreliable" |
| "Ks < 0" | Alignment error | Check codon alignment |
