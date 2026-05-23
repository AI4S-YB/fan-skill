# DESeq2

**Goal:** Differential expression analysis from count data using negative binomial GLM
**Best for:** Standard RNA-seq with ≥3 biological replicates per group

## Prerequisites
- R 4.0+, DESeq2 (Bioconductor)
- Raw count matrix (genes × samples)
- Sample metadata table (condition/group assignments)

## Basic Usage

```r
library(DESeq2)

# Prepare data
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition
)

# Pre-filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results for a contrast
res <- results(dds, contrast = c("condition", "treatment", "control"),
               alpha = 0.05)
res <- res[order(res$padj), ]
summary(res)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| design | ~ condition | Simple group comparison |
| design | ~ batch + condition | Batch correction |
| alpha | 0.05 | FDR threshold |
| lfcThreshold | 0 (default) | Set to 1 for biological significance filtering |

## Plant-Specific Notes

- Plant RNA-seq often has more zeros than mammalian — DESeq2's shrinkage handles this well
- For polyploid species, counts include homeologs — same gene name may appear multiple times (one per subgenome)
- Use `tximport` for Salmon/Kallisto transcript-level quantification

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "model matrix not full rank" | Too many covariates, not enough replicates | Simplify design formula |
| "all counts are zero" | Wrong gene annotation | Check GFF/GTF file used for quantification |
| "contrast not found" | Condition name mismatch | Check `levels(dds$condition)` |
