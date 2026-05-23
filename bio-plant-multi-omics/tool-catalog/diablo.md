# DIABLO — Data Integration Analysis for Biomarker discovery using Latent cOmponents

**Goal:** Supervised multi-omics integration that maximizes class separation while identifying key discriminative features.

**Best for:** Presence of class labels (treatment groups, phenotype categories, varieties), classification with biological interpretability.

**R package:** mixOmics (Bioconductor)

## Prerequisites

- R 4.0+, Bioconductor 3.12+
- Packages: mixOmics, data.table, ggplot2, igraph

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")
```

## Basic Usage

```r
library(mixOmics)
library(data.table)

# ---- Step 1: Load and prepare data ----
# Each omics is a matrix: samples (rows) x features (columns)
# Class labels must be a factor

rna    <- as.matrix(fread("transcriptome_samples_x_genes.csv"), rownames = 1)
met    <- as.matrix(fread("metabolome_samples_x_metabolites.csv"), rownames = 1)
labels <- factor(fread("sample_labels.csv")$condition)

# Verify sample order matches
stopifnot(all(rownames(rna) == rownames(met)))
stopifnot(all(rownames(rna) == names(labels)))

# Combine omics into a list
X <- list(
  Transcriptome = rna,
  Metabolome   = met
)
Y <- labels

# ---- Step 2: Design matrix ----
# Defines relationships between omics blocks
# 0 = no relationship, 0.1 = weak, 0.5 = moderate, 1.0 = full
# For 2 omics, the design is a 2x2 matrix
n_omics <- length(X)
design <- matrix(0.1, ncol = n_omics, nrow = n_omics,
                 dimnames = list(names(X), names(X)))
diag(design) <- 0  # No self-relationships

# ---- Step 3: Tune keepX (number of features to retain per component) ----
set.seed(42)
tune_res <- tune.block.splsda(
  X = X, Y = Y,
  ncomp = 3,
  test.keepX = list(
    Transcriptome = c(5, 10, 20, 50, 100),
    Metabolome    = c(5, 10, 20, 50)
  ),
  design = design,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  dist = "max.dist"
)

# Best keepX values
list.keepX <- tune_res$choice.keepX

# ---- Step 4: Run DIABLO ----
ncomp <- 3  # Number of components
diablo_res <- block.splsda(
  X = X, Y = Y,
  ncomp = ncomp,
  keepX = list.keepX,
  design = design
)

# ---- Step 5: Evaluate performance ----
# Classification error rate
perf_diablo <- perf(diablo_res, validation = "Mfold",
                    folds = 5, nrepeat = 10, dist = "max.dist")
plot(perf_diablo)

# Select optimal number of components
perf_diablo$error.rate

# ---- Step 6: Visualize results ----
# Sample plot (first two components)
plotIndiv(diablo_res, comp = 1:2,
          ind.names = FALSE, legend = TRUE,
          title = "DIABLO Sample Plot")

# Variable contribution plot
plotVar(diablo_res, var.names = c(FALSE, FALSE),
        legend = TRUE, title = "DIABLO Correlation Circle")

# Circos plot — cross-omics feature correlations
circosPlot(diablo_res, cutoff = 0.7,
           line = TRUE, color.blocks = c("darkorchid", "brown1"),
           color.cor = c("chocolate", "grey20"),
           title = "Cross-Omics Correlations")

# Loading plot
plotLoadings(diablo_res, comp = 1,
             contrib = "max", method = "median")

# ---- Step 7: Export results ----
# Selected features per component
selected_vars <- selectVar(diablo_res, comp = 1)
write.csv(selected_vars$Transcriptome$value, "outputs/DIABLO_comp1_RNA_features.csv")
write.csv(selected_vars$Metabolome$value, "outputs/DIABLO_comp1_Met_features.csv")

# Sample scores
write.csv(diablo_res$variates$Transcriptome, "outputs/DIABLO_sample_scores.csv")
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|-------------|-----------|
| design (off-diagonal) | 0.1 | Weak association prevents overfitting; higher values (0.5+) risk overfitting to noise |
| ncomp | 2-5 | More components risk overfitting; use perf() error rate to decide |
| keepX | tuned via tune() | Automatic selection via repeated CV; do not guess |
| validation | "Mfold" | M-fold CV is robust for moderate sample sizes |
| nrepeat | 10-50 | More repeats = more stable error estimate (50 for final, 10 for tuning) |
| dist | "max.dist" | Standard distance; "centroids.dist" for very imbalanced classes |

## Plant-Specific Notes

- **Polyploid species**: Class labels may be subgenome origin or ploidy level. DIABLO can identify subgenome-biased features — check whether selected features are enriched in a particular subgenome.
- **Plant treatment groups**: Drought, salt, nutrient deficiency, pathogen — DIABLO excels at finding the molecular signature of stress responses across omics.
- **Multi-environment**: If samples come from different environments, add environment as a DIABLO outcome OR use it to define categorical labels (environment-specific response).
- **Non-model species**: Use homology-based annotation for biological interpretation. Selected features without annotation may still be valid — flag for follow-up.
- **Developmental stages**: Time-series multi-omics with stage labels — DIABLO can rank features by their contribution to distinguishing developmental stages.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Error in check_Data" | Sample order mismatch | Verify `rownames(X[[1]])` matches all blocks |
| "design matrix not square" | Wrong design dimensions | Must be n_omics x n_omics |
| "too many NAs" | Missing values in input | Impute missing values before DIABLO |
| tune.block.splsda takes forever | Too many keepX values tested | Limit test.keepX to 4-5 values per omics |
| All components have same error | ncomp too high | Reduce ncomp; use perf() to find the elbow |
