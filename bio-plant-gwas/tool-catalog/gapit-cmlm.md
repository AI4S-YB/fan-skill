# GAPIT CMLM (Compressed Mixed Linear Model)

**Goal:** GWAS with compressed MLM — controlling both population structure and kinship
**Best for:** Inbred species, small populations (<200), low-to-medium density markers

## Prerequisites
- R 4.0+, GAPIT 3.0+
- Genotype: Hapmap or numerical format
- Phenotype: CSV/TXT, first column = sample ID

## Basic Usage

```r
library(GAPIT3)

myGAPIT <- GAPIT(
  G = myGD,              # Genotype data (Hapmap format)
  GD = myGD,
  GM = myGM,             # SNP map (chr, snp_name, position)
  Y = myY,               # Phenotype (sample_id, trait_value)
  model = "CMLM",
  PCA.total = 3,         # Number of PCs to include
  kinship.cluster = "average",
  group.from = 200,      # Compression: start at 200 groups
  group.to = 1000,       # Compression: end at 1000 groups
  group.by = 100,        # Compression: step size
  cutOff = 0.05          # P-value threshold for output
)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| PCA.total | 3-5 | Fewer = residual stratification, more = over-correction |
| group.from | max(100, n_samples/2) | Start compression from half the samples |
| group.to | n_samples | End at full sample count |
| kinship.cluster | "average" | Standard hierarchical clustering |

## Plant-Specific Notes

- **Inbred species**: CMLM is the default recommendation. K matrix is essential.
- **Outcross species**: Add PCA.total=5 if population structure is the primary concern.
- **Low-density markers**: CMLM handles low density better than FarmCPU/BLINK.
- **Hapmap format**: Numeric genotype coding (0/1/2 for major allele count).

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "GAPIT not found" | Package not installed | `devtools::install_github("jiabowang/GAPIT3")` |
| "subscript out of bounds" | Mismatch between genotype and phenotype samples | Check sample IDs match |
| "NA/NaN in phenotype" | Missing values not handled | Remove or impute before analysis |
