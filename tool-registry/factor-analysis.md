# Multi-Omics Factor Interpretation

**Goal:** Interpret latent factors from MOFA2 — variance decomposition, top feature extraction, enrichment analysis, and biological annotation.

**Best for:** Post-MOFA2 analysis; understanding what each factor represents biologically.

**R packages:** MOFA2, clusterProfiler, topGO, data.table, ggplot2

## Prerequisites

- R 4.0+
- Packages: MOFA2, clusterProfiler (Bioconductor), topGO (Bioconductor), data.table, ggplot2, pheatmap

## Factor Interpretation Pipeline

```r
library(MOFA2)
library(data.table)
library(ggplot2)
library(pheatmap)

# ---- Load trained MOFA model ----
model <- load_model("outputs/MOFA2_model.hdf5")

# ============================================================
# Step 1: Overall variance decomposition
# ============================================================

# Variance explained by each factor, stacked by view (omics)
plot_variance_explained(model, x = "factor", y = "variance_explained",
                        plot_total = TRUE) +
  ggtitle("Variance Explained per Factor (Stacked by Omics)")

# Total variance explained (all factors combined) — per view
plot_variance_explained(model, x = "view", y = "variance_explained",
                        plot_total = FALSE) +
  ggtitle("Total Variance Explained per Omics Layer")

# Extract numeric values
var_explained <- get_variance_explained(model)
print(var_explained$r2_per_factor)    # R2 per factor per view
print(var_explained$r2_total)         # Total R2 per view

# Check cumulative variance of top 5 factors
top5_r2 <- sum(head(var_explained$r2_total[[1]], 5))
cat(sprintf("Cumulative variance explained by top 5 factors: %.1f%%\n",
            100 * top5_r2))

# ============================================================
# Step 2: Extract top features per factor per view
# ============================================================

# Get weights matrix
weights <- get_weights(model)

# Function: extract top N features (by absolute weight) per factor per view
get_top_features <- function(weights, factor_idx, view_name, n = 20) {
  w <- weights[[view_name]][, factor_idx]
  top_idx <- order(abs(w), decreasing = TRUE)[1:min(n, length(w))]
  data.frame(
    feature = names(w)[top_idx],
    weight = w[top_idx],
    abs_weight = abs(w[top_idx]),
    stringsAsFactors = FALSE
  )
}

# Extract top 20 features for Factor 1 from each omics layer
for (view_name in views_names(model)) {
  cat(sprintf("\n=== Factor 1, View: %s ===\n", view_name))
  top_feat <- get_top_features(weights, factor_idx = 1, view_name = view_name, n = 20)
  print(top_feat)
  write.csv(top_feat,
            sprintf("outputs/factor1_%s_top20_features.csv", view_name),
            row.names = FALSE)
}

# ============================================================
# Step 3: Heatmap of top feature weights across factors
# ============================================================

# For a specific view, show top 30 features (by max absolute weight) across all factors
plot_top_weights_across_factors <- function(model, view_name, n_features = 30, n_factors = 5) {
  w <- get_weights(model)[[view_name]]

  # Select top features by max absolute weight across factors
  max_abs_weight <- apply(w[, 1:min(n_factors, ncol(w)), drop = FALSE], 1,
                          function(x) max(abs(x)))
  top_features <- names(sort(max_abs_weight, decreasing = TRUE)[1:n_features])

  # Subset weight matrix
  w_sub <- w[top_features, 1:min(n_factors, ncol(w)), drop = FALSE]

  # Plot heatmap
  pheatmap(w_sub,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           main = sprintf("Top %d Features — %s (Factors 1-%d)",
                          n_features, view_name, n_factors),
           fontsize_row = 6,
           filename = sprintf("outputs/heatmap_%s_top%d_factors.pdf",
                              view_name, n_factors))
}

for (view_name in views_names(model)) {
  plot_top_features_across_factors(model, view_name, n_features = 30, n_factors = 5)
}

# ============================================================
# Step 4: Factor correlation with sample metadata
# ============================================================

# If you have sample metadata (treatment, tissue, genotype, etc.)
# Correlate factor values with metadata to interpret factors

factor_scores <- get_factors(model)[[1]]  # samples x factors

# Example: read metadata
metadata <- fread("sample_metadata.csv", data.table = FALSE)
rownames(metadata) <- metadata$sample_id
metadata <- metadata[rownames(factor_scores), ]

# Correlate each factor with continuous metadata variables
factor_correlations <- cor(factor_scores,
                           metadata[, sapply(metadata, is.numeric)],
                           use = "pairwise.complete.obs")
print(factor_correlations)

# ANOVA: factor ~ categorical metadata
for (meta_col in colnames(metadata)) {
  if (is.factor(metadata[[meta_col]]) || is.character(metadata[[meta_col]])) {
    cat(sprintf("\n=== Factor ~ %s ===\n", meta_col))
    for (f in 1:ncol(factor_scores)) {
      aov_res <- summary(aov(factor_scores[, f] ~ as.factor(metadata[[meta_col]])))
      p_val <- aov_res[[1]]$`Pr(>F)`[1]
      if (p_val < 0.05) {
        cat(sprintf("  Factor %d: p = %.4f *\n", f, p_val))
      }
    }
  }
}

# ============================================================
# Step 5: GO/KEGG enrichment for top factor features
# ============================================================

library(clusterProfiler)

# Example: GO enrichment for transcriptome top features in Factor 1
run_go_enrichment <- function(gene_list, background_genes, org_db = "org.At.tair.db") {
  library(org_db, character.only = TRUE)

  ego <- enrichGO(
    gene          = gene_list,
    universe      = background_genes,
    OrgDb         = get(org_db),
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  if (nrow(as.data.frame(ego)) > 0) {
    dotplot(ego, showCategory = 20,
            title = "GO Biological Process Enrichment")
  }
  ego
}

# ---- Example usage (requires species-specific OrgDb) ----
# top_genes_factor1 <- get_top_features(weights, 1, "Transcriptome", 100)$feature
# all_genes <- rownames(weights[["Transcriptome"]])
# go_res <- run_go_enrichment(top_genes_factor1, all_genes, "org.At.tair.db")

# ============================================================
# Step 6: Factor interpretation summary
# ============================================================

# Generate a summary table for each factor
summarize_factors <- function(model, n_factors = 5) {
  var_exp <- get_variance_explained(model)$r2_per_factor
  weights <- get_weights(model)

  for (f in 1:n_factors) {
    cat(sprintf("\n========== Factor %d ==========\n", f))

    for (view_name in views_names(model)) {
      r2 <- var_exp[[view_name]][f]
      cat(sprintf("  %s: R2 = %.3f\n", view_name, r2))

      # Top 5 positive and negative weights
      w <- weights[[view_name]][, f]
      top_pos <- names(sort(w, decreasing = TRUE)[1:5])
      top_neg <- names(sort(w, decreasing = FALSE)[1:5])
      cat(sprintf("    Top positive: %s\n", paste(top_pos, collapse = ", ")))
      cat(sprintf("    Top negative: %s\n", paste(top_neg, collapse = ", ")))
    }
  }
}

summarize_factors(model, n_factors = 5)
```

## Plant-Specific Interpretation Notes

- **Factor enriched in "unknown" metabolites**: Do not dismiss. These factors may represent uncharacterized specialized metabolic pathways. Cross-reference with retention time and MS/MS fragmentation patterns.

- **Factor driven by chloroplast genes**: Almost always represents photosynthetic activity variation across samples. Check correlation with sampling time and tissue type.

- **Factor dominated by stress-responsive TFs**: Common in plant datasets. MYB, WRKY, NAC, bHLH families are heavily represented. Check if factor scores correlate with treatment conditions.

- **Polyploid interpretation**: If a factor loads strongly on homeolog pairs (e.g., TraesCS1A... and TraesCS1B... both with high weight), this suggests subgenome-coordinated regulation. If only one homeolog loads, this suggests subgenome-dominant expression.

- **Plant-specific GO terms**: Use `clusterProfiler` with appropriate OrgDb. For non-model species, use `enricher()` with custom gene-to-GO mapping built from InterProScan/eggNOG-mapper results.
