# Multi-Omics Data Preprocessing

**Goal:** Prepare multi-omics data for integration — filter to common samples, normalize per-omics, handle missing values, and batch correction.

**Best for:** All multi-omics integration methods (MOFA2, DIABLO, mixOmics)

**R packages:** data.table, limma, sva, DESeq2, vsn

## Prerequisites

- R 4.0+
- Packages: data.table, limma, sva, DESeq2, vsn, preprocessCore

## Complete Preprocessing Pipeline

```r
library(data.table)

# ============================================================
# Step 1: Load all omics data
# ============================================================
# Each omics should be a matrix: features (rows) x samples (columns)
# For DIABLO/mixOmics, transpose later to samples x features

rna  <- as.matrix(fread("transcriptome_raw_counts.csv"), rownames = 1)
met  <- as.matrix(fread("metabolome_intensities.csv"), rownames = 1)
prot <- as.matrix(fread("proteome_lfq.csv"), rownames = 1)

omics_list <- list(
  Transcriptome = rna,
  Metabolome   = met,
  Proteome     = prot
)

# ============================================================
# Step 2: Filter to common samples
# ============================================================
# Keep only samples present in ALL omics layers

common_samples <- Reduce(intersect, lapply(omics_list, colnames))
cat(sprintf("Common samples across omics: %d\n", length(common_samples)))

if (length(common_samples) < 10) {
  warning("Very few common samples — integration may be unreliable")
}

# Subset each omics to common samples
omics_filtered <- lapply(omics_list, function(m) m[, common_samples, drop = FALSE])

# ============================================================
# Step 3: Feature filtering (per omics)
# ============================================================

# --- Transcriptome: Filter lowly expressed genes ---
filter_transcriptome <- function(count_matrix) {
  # Keep genes with count >= 10 in at least 20% of samples
  keep <- rowSums(count_matrix >= 10) >= (0.2 * ncol(count_matrix))
  cat(sprintf("Transcriptome: %d / %d genes retained after filtering\n",
              sum(keep), nrow(count_matrix)))
  count_matrix[keep, , drop = FALSE]
}

# --- Metabolome: Filter features with too many missing values ---
filter_metabolome <- function(met_matrix, max_missing_frac = 0.5) {
  missing_frac <- rowSums(is.na(met_matrix)) / ncol(met_matrix)
  keep <- missing_frac <= max_missing_frac
  cat(sprintf("Metabolome: %d / %d features retained after filtering\n",
              sum(keep), nrow(met_matrix)))
  met_matrix[keep, , drop = FALSE]
}

# --- Proteome: Filter proteins with too many missing values ---
filter_proteome <- function(prot_matrix, max_missing_frac = 0.3) {
  missing_frac <- rowSums(is.na(prot_matrix) | prot_matrix == 0) / ncol(prot_matrix)
  keep <- missing_frac <= max_missing_frac
  cat(sprintf("Proteome: %d / %d proteins retained after filtering\n",
              sum(keep), nrow(prot_matrix)))
  prot_matrix[keep, , drop = FALSE]
}

omics_filtered$Transcriptome <- filter_transcriptome(omics_filtered$Transcriptome)
omics_filtered$Metabolome   <- filter_metabolome(omics_filtered$Metabolome)
omics_filtered$Proteome     <- filter_proteome(omics_filtered$Proteome)

# ============================================================
# Step 4: Impute missing values (per omics)
# ============================================================

# --- Metabolome: Impute with half-minimum ---
impute_halfmin <- function(mat) {
  mat_imputed <- mat
  for (i in 1:nrow(mat)) {
    row_vals <- mat[i, ]
    if (any(is.na(row_vals))) {
      half_min <- min(row_vals, na.rm = TRUE) / 2
      mat_imputed[i, is.na(row_vals)] <- half_min
    }
  }
  mat_imputed
}

# --- Proteome: Impute with minimum per protein ---
impute_min <- function(mat) {
  mat_imputed <- mat
  for (i in 1:nrow(mat)) {
    row_vals <- mat[i, ]
    nas <- is.na(row_vals) | row_vals == 0
    if (any(nas)) {
      mat_imputed[i, nas] <- min(row_vals[!nas], na.rm = TRUE)
    }
  }
  mat_imputed
}

omics_filtered$Metabolome <- impute_halfmin(omics_filtered$Metabolome)
omics_filtered$Proteome   <- impute_min(omics_filtered$Proteome)

# ============================================================
# Step 5: Normalize (per omics)
# ============================================================

# --- Transcriptome: VST (DESeq2) or log2(CPM+1) ---
normalize_rna_vst <- function(count_matrix) {
  library(DESeq2)
  # Create minimal DESeq2 object
  col_data <- data.frame(
    row.names = colnames(count_matrix),
    condition = factor(rep("A", ncol(count_matrix)))
  )
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = col_data,
    design = ~ 1
  )
  # Variance stabilizing transformation
  vsd <- vst(dds, blind = TRUE)
  assay(vsd)
}

# --- Metabolome: Log2 + Pareto scaling ---
normalize_metabolome <- function(met_matrix) {
  # Log2 transform
  mat_log <- log2(met_matrix + 1)

  # Pareto scaling: (x - mean) / sqrt(sd)
  mat_pareto <- t(scale(t(mat_log), center = TRUE, scale = TRUE))
  mat_pareto <- mat_pareto * sqrt(apply(mat_log, 1, sd, na.rm = TRUE))
  # Avoid: Pareto scaling in practice is mean-center + divide by sqrt(sd)
  # Simpler implementation:
  row_means <- rowMeans(mat_log, na.rm = TRUE)
  row_sds   <- apply(mat_log, 1, sd, na.rm = TRUE)
  mat_pareto2 <- (mat_log - row_means) / sqrt(row_sds)
  mat_pareto2
}

# --- Proteome: Quantile normalization ---
normalize_proteome <- function(prot_matrix) {
  library(preprocessCore)
  normalize.quantiles(prot_matrix, keep.names = TRUE)
}

omics_norm <- omics_filtered
omics_norm$Transcriptome <- normalize_rna_vst(omics_filtered$Transcriptome)
omics_norm$Metabolome     <- normalize_metabolome(omics_filtered$Metabolome)
omics_norm$Proteome       <- normalize_proteome(omics_filtered$Proteome)

# ============================================================
# Step 6: Batch correction (optional — only if batch effect detected)
# ============================================================

# Check for batch effect: PCA colored by batch
check_batch_effect <- function(expr_matrix, batch_labels) {
  # Simple PCA
  pca_res <- prcomp(t(expr_matrix), center = TRUE, scale. = TRUE)

  # Plot PC1 vs PC2 colored by batch
  plot(pca_res$x[, 1], pca_res$x[, 2],
       col = as.factor(batch_labels),
       pch = 19, cex = 1.5,
       xlab = "PC1", ylab = "PC2",
       main = "PCA Colored by Batch")

  # ANOVA of PC1 by batch
  aov_pc1 <- summary(aov(pca_res$x[, 1] ~ as.factor(batch_labels)))
  cat(sprintf("PC1 variance explained by batch: p = %.4f\n",
              aov_pc1[[1]]$`Pr(>F)`[1]))
}

# Apply ComBat if needed
apply_combat <- function(expr_matrix, batch_labels) {
  library(sva)
  mod <- model.matrix(~ 1, data = data.frame(row.names = colnames(expr_matrix)))
  combat_res <- ComBat(dat = expr_matrix, batch = batch_labels, mod = mod)
  combat_res
}

# ---- Example usage (commented out) ----
# batch_labels <- fread("sample_batch.csv")$batch
# check_batch_effect(omics_norm$Transcriptome, batch_labels)
# If batch effect significant, run:
# omics_norm$Transcriptome <- apply_combat(omics_norm$Transcriptome, batch_labels)

# ============================================================
# Step 7: Save processed data
# ============================================================

for (name in names(omics_norm)) {
  write.csv(omics_norm[[name]],
            file = sprintf("outputs/%s_preprocessed.csv", name))
  cat(sprintf("Saved %s: %d features x %d samples\n",
              name, nrow(omics_norm[[name]]), ncol(omics_norm[[name]])))
}

cat("\nPreprocessing complete. Ready for integration.\n")
```

## Plant-Specific Preprocessing Notes

- **Chloroplast/Mito transcripts**: PCR duplicates from poly-A selection do not remove organellar RNA. These can dominate count matrices and inflate variance. Check the top 100 most expressed "genes" — if dominated by chloroplast transcripts, consider removing them.

- **Metabolome annotation level**: Plant metabolomes have many "unknown" features (often > 50%). Do NOT filter them out — they carry biological information. Label them as "Unknown_mz_rt" (mass/retention-time).

- **Secondary metabolites**: Many plant specialized metabolites follow bimodal distributions (present vs absent across genotypes). Log transformation with +1 pseudocount handles this, but check the distribution before and after.

- **Seasonal/tissue variation**: Multi-tissue plant datasets have inherently different baseline expression profiles. Standard normalization per omics is sufficient; do not quantile-normalize across tissues.

- **Polyploid mapping bias**: RNA-seq read mapping bias toward one subgenome can create spurious signals. If suspect, check expression levels of known single-copy genes across subgenomes.

## QC Output Table

| Checkpoint | What to Check | Pass Threshold |
|------------|---------------|----------------|
| Common samples | Number retained | >= 50% of minimum per omics |
| Feature filtering | Features per omics | >= 100 |
| Missing values | % imputed per feature | < 50% |
| PCA post-processing | PC1 by batch | p > 0.05 |
| Expression distribution | Density plots | Similar across samples |
