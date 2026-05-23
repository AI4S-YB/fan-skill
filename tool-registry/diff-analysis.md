# Differential Metabolite Analysis (limma)

**Goal:** Identify metabolites significantly different between experimental conditions
**Best for:** Metabolomics studies with >= 3 biological replicates per group

## Prerequisites
- R 4.0+ with limma, ggplot2, pheatmap packages
- Processed feature table (log2-transformed peak areas)
- Sample metadata with group assignments

## limma Workflow

### 1. Load and Prepare Data

```r
library(limma)

# Read feature table and metadata
feature_table <- read.csv("feature_table_log2.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1)

# Verify sample order matches
stopifnot(all(colnames(feature_table) == rownames(metadata)))

# Define experimental design
group <- factor(metadata$group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# For paired design (e.g., same plant before/after treatment)
# design <- model.matrix(~ 0 + group + donor)
```

### 2. Fit Linear Model

```r
# Fit model
fit <- lmFit(feature_table, design)

# Define contrasts
contrast.matrix <- makeContrasts(
  Treatment_vs_Control = Treatment - Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)  # trend=TRUE recommended for metabolomics
```

### 3. Extract Results

```r
# Get all results
results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH",
                     sort.by = "p")

# Add feature metadata (m/z, RT, putative annotation)
feature_info <- read.csv("feature_info.csv", row.names = 1)
results <- merge(results, feature_info, by = "row.names", all.x = TRUE)

# Flag significant features
results$significant <- results$adj.P.Val < 0.05
results$direction <- ifelse(results$logFC > 0, "Up", "Down")

# Export
write.csv(results, "differential_metabolites.csv", row.names = FALSE)

# Significant subset
sig_results <- subset(results, adj.P.Val < 0.05)
write.csv(sig_results, "significant_metabolites.csv", row.names = FALSE)

cat(sprintf("Total features: %d\n", nrow(results)))
cat(sprintf("Significant (q < 0.05): %d\n", nrow(sig_results)))
cat(sprintf("  Up: %d, Down: %d\n",
            sum(sig_results$direction == "Up"),
            sum(sig_results$direction == "Down")))
```

### 4. QC Visualizations

```r
# Volcano plot
library(ggplot2)
ggplot(results, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("grey50", "#C73E1D")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot: Treatment vs Control") +
  theme_bw(base_size = 14)

# Mean-difference (MA) plot
plotMD(fit2, column = 1, status = results$significant,
       main = "MA Plot: Treatment vs Control")
abline(h = 0, col = "red", lty = 2)
```

### 5. Heatmap of Significant Metabolites

```r
library(pheatmap)

# Extract significant feature intensities
sig_matrix <- feature_table[rownames(feature_table) %in% sig_results$Row.names, ]

# Z-score normalize rows
sig_matrix_z <- t(scale(t(sig_matrix)))

# Annotation
annotation_col <- data.frame(
  Group = group,
  row.names = colnames(sig_matrix_z)
)

pheatmap(sig_matrix_z,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("#2E86AB", "white", "#C73E1D"))(100),
         main = "Significant Metabolites (q < 0.05)")
```

## Missing Value Handling

```r
# Before limma: handle missing values
# Option 1: Replace zeros with small value (for left-censored missing)
feature_table[feature_table == 0] <- NA
feature_table <- apply(feature_table, 2, function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE) / 2
  return(x)
})

# Option 2: Impute using k-nearest neighbors
library(impute)
feature_table_imputed <- impute.knn(as.matrix(feature_table))$data

# Option 3: Filter features with too many missing values
max_missing <- 0.5  # Keep features present in >50% of samples
keep <- rowSums(is.na(feature_table)) / ncol(feature_table) < max_missing
feature_table_filtered <- feature_table[keep, ]
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| trend | Account for mean-variance relationship (TRUE for metabolomics) |
| adjust.method | Multiple testing correction (BH recommended) |
| robust | Robust regression for outliers |
| logFC cutoff | Biological significance threshold (>= 1 for LC-MS, >= 0.58 for GC-MS) |

## Plant-Specific Considerations

- Plant metabolites often show high biological variability — 4-6 replicates recommended
- Diurnal variation: time-of-day sampling must be standardized; include time as covariate
- Developmental stage: young vs mature leaves produce different metabolomes; match tissue age
- Stress treatments: wounding during sampling triggers rapid metabolic changes (seconds to minutes); standardize collection protocol

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| No significant features | Small effect size or low power | Increase MC samples; lower q threshold to 0.1 |
| "Design matrix not full rank" | Confounded design variables | Simplify formula; check metadata |
| Too many significant features | Batch effect confounded with condition | Check experimental design; use batch correction |
| All logFC near zero | Data not properly normalized | Use QC-RLSC or ComBat before analysis |
