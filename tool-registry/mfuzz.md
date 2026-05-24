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

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| m (fuzzifier) | 1.5 | For sharp on/off patterns (stress response), lower to 1.1-1.3; for graded developmental transitions, raise to 1.8-2.0 | Controls cluster boundary softness; low m approaches hard clustering (cleaner patterns), high m allows genes to bridge clusters (captures biological gradients) |
| centers (cluster count) | `mestimate()` auto-estimate | If mestimate returns 1 (all genes similar), manually set 4-8 and pre-filter to top 50% variable genes | Auto-estimation fails when expression patterns are homogeneous; filtering to variable genes reveals hidden sub-patterns |
| membership threshold | 0.5 | For core regulon definition, raise to 0.7; for exploratory GO enrichment per cluster, lower to 0.3 | Higher thresholds yield smaller, highly confident gene sets; lower thresholds include transitional genes that may have dual functions |
| standardize method | z-score (gene-wise) | For absolute expression level comparisons (e.g., cross-species), consider using log2-transformed raw values without standardization | Standardization removes baseline expression differences; preserved when comparing absolute levels across conditions |
| Minimum time points | >=5 | With 3-4 time points, reduce m to 1.1-1.2 and use fewer centers (k=3-6) | Mfuzz cluster shapes become ambiguous with few time points; soft clustering with limited temporal resolution can produce unstable memberships |

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
- **Time point count sensitivity**: Mfuzz cluster reliability drops with fewer than 5 time points. The fuzzifier parameter (m) cannot compensate for sparse temporal sampling — with only 3 time points, even m=2.0 cannot distinguish linear from transient patterns. For dense time courses (>12 points), Mfuzz may over-cluster into too many subtly different groups — increase m to 1.8-2.0 and visually merge similar cluster centers.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| `mestimate returns 1` | All genes have similar pattern | Filter to most variable genes (top 50% by variance) |
| Cluster centers nearly identical | m too high (overly soft) | Reduce m to 1.2; re-run |
| `mfuzz: NA in exprs` | Unfiltered genes with NA | Remove rows with any NA |
| Membership scores uniformly low (<0.3) | Too many clusters requested | Reduce centers or increase m |
| `standardise` creates Inf values | Zero-variance genes present | Filter genes with sd = 0 before standardise() |
