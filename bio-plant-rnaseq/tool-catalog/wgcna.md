# WGCNA (Weighted Gene Co-expression Network Analysis)

**Goal:** Identify co-expressed gene modules and link them to traits
**Best for:** >15 samples with continuous trait measurements

## Prerequisites
- R 4.0+, WGCNA
- Normalized expression data (VST or rlog from DESeq2, or TPM matrix)
- Trait data (optional, for module-trait association)

## Basic Usage

```r
library(WGCNA)
allowWGCNAThreads()

# Choose soft threshold
powers <- c(1:20)
sft <- pickSoftThreshold(t(datExpr), powerVector = powers)
# Select power where scale-free R² > 0.80

# Network construction
net <- blockwiseModules(
  datExpr = t(datExpr),
  power = sft$powerEstimate,
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = FALSE
)

# Module-trait association
moduleTraitCor <- cor(net$MEs, trait_data)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| minModuleSize | 30 (15-60 samples) | Smaller for fewer genes |
| mergeCutHeight | 0.25 | Module merging threshold |
| power | auto (pickSoftThreshold) | Based on scale-free topology fit |

## Plant-Specific Notes

- Plant co-expression modules often reflect tissue/organ-specific patterns
- Modules correlated with yield/quality traits are candidates for functional validation
- Polyploid homeologs may co-express in the same module (if retained together) or separate modules (if subfunctionalized)
