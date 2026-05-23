# Mfuzz: Soft Clustering for Time-Series Expression

**Goal:** Soft-cluster genes by temporal expression pattern, allowing each
gene partial membership in multiple clusters.

**Best for:** Large gene sets (>= 1000 genes) with complex, graded
expression patterns typical of plant transcriptomes.

## Prerequisites

- R 4.0+, Bioconductor `Mfuzz` package
- Smoothed, normalized expression matrix (gene x time-point)
- Gene filtering: remove genes with < 20% of max expression or zero variance

## Basic Usage

```r
library(Mfuzz)

# 1. Create ExpressionSet
eset <- new("ExpressionSet", exprs = as.matrix(smoothed_mat))

# 2. Standardize (gene-wise z-score)
eset <- standardise(eset)

# 3. Estimate optimal cluster number
m1 <- mestimate(eset)
cat("Estimated cluster number:", m1, "\n")

# 4. Alternative: manual selection via cselection
set.seed(42)
cselection_result <- cselection(
  eset,
  m = 1:as.integer(m1 * 2),
  crange = seq(2, 20, 2)
)

# 5. Run Mfuzz clustering
set.seed(42)
cl <- mfuzz(eset, centers = m1, m = 1.5)
# m = fuzzification parameter: 1.0 = hard, > 1.5 = very soft

# 6. Extract membership scores
membership <- cl$membership  # gene x cluster matrix

# 7. Assign genes to primary cluster (highest membership)
primary_cluster <- apply(membership, 1, which.max)
primary_score <- apply(membership, 1, max)

# 8. Filter for confident membership
confident_genes <- names(primary_score[primary_score > 0.5])
cat(length(confident_genes), "genes with membership > 0.5\n")

# 9. Plot cluster centers
pdf("mfuzz_clusters.pdf", width = 10, height = 8)
mfuzz.plot2(
  eset,
  cl = cl,
  mfrow = c(ceiling(m1/3), 3),
  time.labels = colnames(smoothed_mat),
  xlab = "Time point",
  ylab = "Standardized expression"
)
dev.off()
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| centers | `mestimate()` | Number of clusters; auto-estimated or manual |
| m | 1.5 | Fuzzification parameter; 1.0-2.0, higher = softer |
| mfuzz version | standardise() | Z-score standardize (recommended) vs raw |
| threshold | 0.5 | Minimum membership score for "confident" assignment |
| seed | -- | Set for reproducibility; Mfuzz is non-deterministic |

## Plant-Specific Notes

- **Paralog expression**: soft clustering captures the graded divergence
  of paralog expression patterns after whole-genome duplication. Hard
  clustering would force paralogs into separate bins.
- **Tissue atlases**: with many organs, m=1.2-1.5 is sufficient.
  For stress time courses with sharp on/off patterns, lower m (1.1-1.3)
  gives crisper clusters.
- **Rice/Arabidopsis**: typical cluster counts are 8-16 for whole-genome
  developmental series. Stress-specific subsets usually cluster into 4-8 groups.
- **Maize**: with subgenome expression divergence, run Mfuzz on each
  subgenome separately, then compare cluster assignments.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| `mestimate returns 1` | All genes have similar pattern | Filter to most variable genes (top 50% by variance) |
| Cluster centers nearly identical | m too high (overly soft) | Reduce m to 1.2; re-run |
| `mfuzz: NA in exprs` | Unfiltered genes with NA | Remove rows with any NA |
| Membership scores uniformly low (<0.3) | Too many clusters requested | Reduce centers or increase m |
| `standardise` creates Inf values | Zero-variance genes present | Filter genes with sd = 0 before standardise() |
