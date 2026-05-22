# ADMIXTURE

**Goal:** Estimate ancestry proportions for each individual
**Best for:** Detecting subpopulations and admixture patterns

## Prerequisites
- ADMIXTURE 1.3+
- LD-pruned PLINK binary files

## Basic Usage

```bash
for K in {1..10}; do
    admixture --cv input_ld.bed $K | tee admixture_${K}.log
done
grep "CV error" admixture_*.log
```

## Key Parameters
| Parameter | Purpose |
|-----------|---------|
| --cv | Cross-validation to select best K |
| -j | Number of threads |

## Plant-Specific Notes
- Inbred crops often show "mixed" ancestry at correct K due to LD, not true admixture
- Cross-validation may favor higher K in structured plant populations

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "Error reading .bed" | PLINK version mismatch | Recreate BED with PLINK 1.9 |
