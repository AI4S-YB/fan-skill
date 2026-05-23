# PLINK PCA

**Goal:** Principal component analysis on genotype data to visualize population structure
**Best for:** All population genetics analyses — first step before ADMIXTURE or GWAS

## Prerequisites
- PLINK 1.9 or 2.0
- LD-pruned genotype data (BED/BIM/FAM)

## Basic Usage

```bash
plink --bfile input_ld --pca 10 --out pca_result
```

Output: `pca_result.eigenvec` (PC coordinates), `pca_result.eigenval` (eigenvalues)

## Key Parameters

| Parameter | Recommended | Purpose |
|-----------|------------|---------|
| --pca | 10 | Number of PCs to compute |
| --bfile | input_ld | LD-pruned input |

## Plant-Specific Notes
- For polyploids, run per-subgenome PCA
- For inbred species, PCs may separate inbred lines cleanly
- Check for outlier samples (>6 SD on any PC)

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "Error: .bed file not found" | Wrong prefix | Verify with `ls prefix.*` |
| All PCs explain <3% | No LD pruning or data is uniform | Check LD pruning step |
