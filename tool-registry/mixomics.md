# mixOmics PLS / sPLS — Pairwise Omics Integration

**Goal:** Integrate exactly two omics layers using (sparse) Partial Least Squares to find latent components maximizing covariance.

**Best for:** Two omics only (e.g., transcriptome + metabolome), identifying correlated feature pairs.

**R package:** mixOmics (Bioconductor)

## Prerequisites

- R 4.0+, Bioconductor 3.12+
- Packages: mixOmics, data.table, ggplot2

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")
```

## Basic Usage — sPLS (Recommended)

```r
library(mixOmics)
library(data.table)

# ---- Step 1: Load two omics matrices ----
# Rows = samples, Columns = features
# Sample order MUST be identical in both matrices

X_rna <- as.matrix(fread("transcriptome_vst.csv"), rownames = 1)    # Predictor
Y_met <- as.matrix(fread("metabolome_log2.csv"), rownames = 1)      # Response

# Verify sample alignment
stopifnot(all(rownames(X_rna) == rownames(Y_met)))

# ---- Step 2: Tune the number of components and features ----
set.seed(42)

# Tune keepX (features to retain in X per component)
tune_x <- tune.spls(X_rna, Y_met,
                    ncomp = 3,
                    test.keepX = c(5, 10, 20, 50, 100, 200),
                    validation = "Mfold",
                    folds = 5,
                    nrepeat = 10,
                    measure = "cor")

# Tune keepY (features to retain in Y per component)
tune_y <- tune.spls(Y_met, X_rna,
                    ncomp = 3,
                    test.keepX = c(5, 10, 20, 50),
                    validation = "Mfold",
                    folds = 5,
                    nrepeat = 10,
                    measure = "cor")

best_keepX <- tune_x$choice.keepX
best_keepY <- tune_y$choice.keepX

# ---- Step 3: Run sPLS ----
ncomp <- min(3, length(best_keepX))
spls_res <- spls(X_rna, Y_met,
                 ncomp = ncomp,
                 keepX = best_keepX,
                 keepY = best_keepY,
                 mode = "regression")

# ---- Step 4: Evaluate performance ----
# Cross-validation correlation
perf_res <- perf(spls_res, validation = "Mfold",
                 folds = 5, nrepeat = 10)
plot(perf_res)

# Correlation between latent components
plot(spls_res$variates$X[, 1], spls_res$variates$Y[, 1],
     xlab = "X component 1", ylab = "Y component 1",
     main = sprintf("r = %.3f",
                    cor(spls_res$variates$X[,1], spls_res$variates$Y[,1])))

# ---- Step 5: Visualize ----
# Sample plot
plotIndiv(spls_res, comp = 1:2,
          ind.names = FALSE,
          title = "sPLS Sample Plot")

# Variable plot — selected features
plotVar(spls_res, comp = 1:2,
        var.names = c(TRUE, TRUE),
        cutoff = 0.5,
        title = "sPLS Feature Correlation")

# Loading plot for X (transcriptome)
plotLoadings(spls_res, comp = 1,
             block = "X", method = "median",
             contrib = "max")

# Loading plot for Y (metabolome)
plotLoadings(spls_res, comp = 1,
             block = "Y", method = "median",
             contrib = "max")

# ---- Step 6: Network of correlated features ----
# Requires igraph
network(spls_res, comp = 1,
        cutoff = 0.6,
        save = "pdf",
        name.save = "outputs/sPLS_network")

# ---- Step 7: Export results ----
# Selected features per component
selected_x <- selectVar(spls_res, comp = 1, block = "X")
selected_y <- selectVar(spls_res, comp = 1, block = "Y")
write.csv(selected_x$value, "outputs/sPLS_comp1_X_features.csv")
write.csv(selected_y$value, "outputs/sPLS_comp1_Y_features.csv")

# Correlation matrix between selected features
write.csv(cim_res$mat, "outputs/sPLS_feature_correlations.csv")
```

## Standard PLS (without sparsity, for quick exploration)

```r
# Quick initial look — no feature selection
pls_res <- pls(X_rna, Y_met,
               ncomp = 5,
               mode = "regression")

# Variance explained
plot(pls_res$explained_variance$X, type = "b",
     xlab = "Component", ylab = "Variance explained in X",
     main = "X variance explained by PLS components")
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|-------------|-----------|
| ncomp | 2-5 | More components capture more covariance; 5 for thorough exploration |
| keepX (transcriptome) | 20-200 | Retain enough genes for biological pathways |
| keepY (metabolome) | 10-100 | Fewer metabolites typically selected |
| mode | "regression" | Asymmetric: X predicts Y (gene expression -> metabolites) |
| mode | "canonical" | Symmetric: bidirectional relationship |
| validation | "Mfold" | Standard for moderate samples; use "loo" for very small n |

## Plant-Specific Notes

- **Gene-to-metabolite mapping**: sPLS is ideal for transcriptome + metabolome integration in plant specialized metabolism studies (flavonoids, alkaloids, terpenoids). The top-loading genes often include known biosynthetic enzymes.
- **Tissue-specific**: Run sPLS per tissue type — leaf, root, seed metabolomes have fundamentally different compositions and correlation structures with the transcriptome.
- **Multi-treatment**: If samples include treatment and control, run sPLS within each group separately and compare selected features (consistent genes across treatments suggest core regulatory mechanisms).
- **Non-model species**: sPLS does not require gene annotation. Selected features can guide annotation prioritization — genes that correlate strongly with known metabolites deserve annotation effort.
- **Polyploids**: Run separately per subgenome if subgenome assignment is available. Compare homeolog selection patterns.
- **Multi-timepoint**: For time-series two-omics data, consider `timeOmics` or run sPLS per timepoint and track how gene-metabolite correlations evolve.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "number of rows must match" | Sample count differs | Filter to common samples |
| "NAs in X" | Missing values | Impute missing values in both matrices |
| No features selected (keepX=0) | tune selected 0 features | Manually set keepX to 10-20 as a starting point |
| tune.spls runs indefinitely | Large dataset + many test.keepX | Reduce test.keepX to 3-4 values |
| cor = 0 between components | No covariance between omics | Check data quality; may need different normalization |
