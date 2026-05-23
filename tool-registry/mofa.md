# MOFA2 — Multi-Omics Factor Analysis

**Goal:** Unsupervised integration of multi-omics data to discover latent factors that capture coordinated variation across omics layers.

**Best for:** Exploratory analysis, dimension reduction, 3+ omics types, no class labels needed.

**R package:** MOFA2 (Bioconductor)

## Prerequisites

- R 4.0+, Bioconductor 3.12+
- Packages: MOFA2, data.table, ggplot2, pheatmap, reticulate (for Python backend)

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MOFA2")
```

## Basic Usage

```r
library(MOFA2)
library(data.table)

# ---- Step 1: Prepare data as a list of matrices ----
# Each omics is a matrix: features (rows) x samples (columns)
# Sample names MUST match across omics layers

rna <- as.matrix(fread("transcriptome_vst.csv"), rownames = 1)
met <- as.matrix(fread("metabolome_log2.csv"), rownames = 1)
prot <- as.matrix(fread("proteome_lfq.csv"), rownames = 1)

# Create the MOFA object
mofa_data <- create_mofa(list(
  "Transcriptome" = rna,
  "Metabolome"   = met,
  "Proteome"     = prot
))

# ---- Step 2: Set data options ----
data_opts <- get_default_data_options(mofa_data)
# For count-based data, consider Poisson likelihood:
# data_opts$likelihoods <- c("Transcriptome" = "gaussian",
#                             "Metabolome"   = "gaussian",
#                             "Proteome"     = "gaussian")

# ---- Step 3: Set model options ----
model_opts <- get_default_model_options(mofa_data)
model_opts$num_factors <- 10   # Start with 10 factors

# ---- Step 4: Set training options ----
train_opts <- get_default_training_options(mofa_data)
train_opts$maxiter <- 5000
train_opts$convergence_mode <- "fast"
train_opts$drop_factor_threshold <- 0.02  # Drop factors explaining < 2%

# ---- Step 5: Prepare and run MOFA ----
mofa_obj <- prepare_mofa(
  object = mofa_data,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Run (this may take minutes to hours depending on data size)
model <- run_mofa(mofa_obj, outfile = "outputs/MOFA2_model.hdf5")

# ---- Step 6: Basic exploration ----
# Variance explained per factor
plot_variance_explained(model, x = "factor", y = "variance_explained")

# Variance explained by omics per factor (stacked bar)
plot_variance_explained(model, x = "factor", y = "variance_explained",
                        plot_total = TRUE)

# Factor values across samples
plot_factor(model, factors = 1:3, color_by = "Factor1")

# Top weights for each factor
plot_top_weights(model, factor = 1, nfeatures = 10)

# ---- Step 7: Export results ----
# Factor scores (samples x factors)
write.csv(get_factors(model)[[1]], "outputs/factor_scores.csv")

# Feature weights (per factor)
write.csv(get_weights(model), "outputs/feature_weights.csv")
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|-------------|-----------|
| num_factors | min(15, n_samples / 3) | Too many = overfit, too few = miss signals |
| maxiter | 5000 | Sufficient for most plant datasets |
| drop_factor_threshold | 0.02 | Factors < 2% variance likely noise |
| likelihood (RNA) | gaussian | For VST/log-normalized counts |
| likelihood (Metabolome) | gaussian | For log-transformed metabolite intensities |
| likelihood (ATAC) | gaussian | For normalized accessibility scores |

## Plant-Specific Notes

- **Polyploid species**: MOFA2 does not natively handle homeologs. Run subgenome-specific MOFA or prefix homeolog pairs with subgenome labels (e.g., "A_TraesCS1A01G000100", "B_TraesCS1B01G000100") and inspect whether they load on the same factor.
- **Chloroplast/mitochondrial genes**: Keep them in the expression matrix but flag them in downstream interpretation — they reflect organellar activity, not nuclear regulation.
- **Unknown metabolites**: Do NOT remove them before MOFA2. Unknown features can form meaningful factors and be retrospectively identified.
- **Multi-tissue**: If samples span tissues (leaf, root, seed), add tissue as a covariate in `model_opts$covariates` or run MOFA per tissue and compare factor structures.
- **Small sample size (n<20)**: Limit `num_factors` to <= 5. Use `train_opts$seed` for reproducibility. Report factor stability via bootstrap.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Sample names do not match" | Different sample IDs across omics | Standardize sample names before creating MOFA object |
| "Python not found" | MOFA2 backend requires Python | `reticulate::py_install("mofapy2", pip = TRUE)` |
| "No factors retained after dropping" | All factors below variance threshold | Lower `drop_factor_threshold` to 0.01 |
| "NaN in factor values" | Zero-variance features | Filter features with variance = 0 before input |
| Memory error | Too many features | Filter to top 5000 variable features per omics |
