# Multi-Omics Visualization

**Goal:** Publication-quality visualization of multi-omics integration results — factor plots, heatmaps, sample ordination, cross-omics networks, and DIABLO diagnostic plots.

**Best for:** All integration results (MOFA2, DIABLO, mixOmics PLS/sPLS)

**R packages:** ggplot2, pheatmap, ComplexHeatmap, igraph, RColorBrewer, cowplot

## Prerequisites

- R 4.0+
- Packages: ggplot2, pheatmap, ComplexHeatmap (Bioconductor), RColorBrewer, cowplot, ggrepel

## Visualization Code

```r
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(ggrepel)

# ============================================================
# 1. Variance Decomposition Plot (MOFA2)
# ============================================================

plot_variance_decomposition <- function(model, outfile = "outputs/figures/variance_decomposition.pdf") {
  # Manual plot for more control than MOFA2 defaults
  var_exp <- get_variance_explained(model)
  r2_per_factor <- var_exp$r2_per_factor

  # Convert to long-format data frame
  df_list <- list()
  for (view_name in names(r2_per_factor)) {
    r2 <- r2_per_factor[[view_name]]
    df_list[[view_name]] <- data.frame(
      Factor = factor(seq_along(r2), levels = seq_along(r2)),
      R2 = r2,
      Omics = view_name
    )
  }
  df <- do.call(rbind, df_list)

  # Calculate total R2 per factor
  total_r2 <- aggregate(R2 ~ Factor, data = df, FUN = sum)

  p <- ggplot(df, aes(x = Factor, y = R2, fill = Omics)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    geom_text(data = total_r2,
              aes(x = Factor, y = R2 + 0.02,
                  label = sprintf("%.1f%%", 100 * R2),
                  fill = NULL),
              size = 3, vjust = 0) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Variance Explained per Factor",
         subtitle = sprintf("Total R2 across all omics: %.1f%%",
                            100 * sum(total_r2$R2)),
         x = "Factor", y = "Variance Explained (R2)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(outfile, p, width = 10, height = 6)
  cat(sprintf("Saved: %s\n", outfile))
}

# ============================================================
# 2. Factor-Feature Heatmap (MOFA2)
# ============================================================

plot_factor_heatmap <- function(model, factor_idx = 1, view_name = NULL,
                                n_features = 50, outfile = NULL) {
  if (is.null(view_name)) {
    view_name <- views_names(model)[1]
  }
  if (is.null(outfile)) {
    outfile <- sprintf("outputs/figures/factor%d_%s_heatmap.pdf",
                       factor_idx, view_name)
  }

  w <- get_weights(model)[[view_name]]
  z <- get_factors(model)[[1]]

  # Select top features for this factor
  top_idx <- order(abs(w[, factor_idx]), decreasing = TRUE)[1:min(n_features, nrow(w))]
  w_top <- w[top_idx, factor_idx, drop = FALSE]

  # Data matrix for selected features (samples x features)
  data_mat <- t(get_data(model, views = view_name)[[view_name]][rownames(w_top), , drop = FALSE])

  # Sort samples by factor value
  sample_order <- order(z[, factor_idx])
  data_mat <- data_mat[sample_order, , drop = FALSE]

  # Scale features for heatmap
  data_mat_scaled <- scale(data_mat)

  # Annotation: factor value
  annotation_col <- data.frame(
    Factor = z[sample_order, factor_idx],
    row.names = rownames(data_mat)
  )

  pheatmap(t(data_mat_scaled),
           annotation_col = annotation_col,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           show_rownames = TRUE,
           fontsize_row = 6,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           main = sprintf("Factor %d Top Features — %s", factor_idx, view_name),
           filename = outfile)

  cat(sprintf("Saved: %s\n", outfile))
}

# ============================================================
# 3. Factor Scores — Sample Ordination (MOFA2)
# ============================================================

plot_factor_scores <- function(model, factor_x = 1, factor_y = 2,
                               color_by = NULL, label_by = NULL,
                               outfile = "outputs/figures/factor_scores.pdf") {
  z <- get_factors(model)[[1]]

  df <- data.frame(
    Sample = rownames(z),
    FactorX = z[, factor_x],
    FactorY = z[, factor_y]
  )

  # Add color groups if provided
  if (!is.null(color_by) && is.character(color_by)) {
    df$Group <- color_by
  } else if (!is.null(color_by) && is.data.frame(color_by)) {
    df$Group <- color_by[rownames(z), 1]
  } else {
    df$Group <- "All"
  }

  p <- ggplot(df, aes(x = FactorX, y = FactorY, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linewidth = 0.5) +
    labs(x = sprintf("Factor %d", factor_x),
         y = sprintf("Factor %d", factor_y),
         title = sprintf("Factor %d vs Factor %d — Sample Scores", factor_x, factor_y)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")

  if (!is.null(label_by)) {
    df$Label <- label_by
    p <- p + geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20)
  }

  ggsave(outfile, p, width = 8, height = 7)
  cat(sprintf("Saved: %s\n", outfile))
}

# ============================================================
# 4. Multi-Omics Correlation Heatmap
# ============================================================

plot_multiomics_correlation <- function(omics_list, outfile = "outputs/figures/omics_correlation.pdf") {
  # Combine selected top features from each omics into a correlation matrix

  # Take top 100 variable features per omics (for computational efficiency)
  top_features_list <- lapply(omics_list, function(mat) {
    vars <- apply(mat, 1, var, na.rm = TRUE)
    names(sort(vars, decreasing = TRUE)[1:min(100, nrow(mat))])
  })

  # Combine all selected features into one matrix
  combined_mat <- do.call(rbind, lapply(names(omics_list), function(name) {
    omics_list[[name]][top_features_list[[name]], , drop = FALSE]
  }))

  # Compute correlation across features
  cor_mat <- cor(t(combined_mat), use = "pairwise.complete.obs")

  # Annotation: which omics each feature belongs to
  annotation_row <- data.frame(
    Omics = factor(rep(names(omics_list), sapply(top_features_list, length))),
    row.names = rownames(cor_mat)
  )

  # Color palette for omics
  omics_colors <- brewer.pal(max(3, length(omics_list)), "Set1")
  names(omics_colors) <- names(omics_list)
  annotation_colors <- list(Omics = omics_colors)

  pheatmap(cor_mat,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           show_colnames = FALSE,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           main = "Cross-Omics Feature Correlation Heatmap",
           filename = outfile,
           width = 12, height = 10)

  cat(sprintf("Saved: %s\n", outfile))
}

# ============================================================
# 5. DIABLO Diagnostic Plots
# ============================================================

plot_diablo_diagnostics <- function(diablo_res, perf_res, out_prefix = "outputs/figures/DIABLO_") {

  # 5a. Classification error rate by component
  pdf(paste0(out_prefix, "error_rate.pdf"), width = 8, height = 6)

  # Plot overall error rate
  overall_err <- perf_res$error.rate$overall
  if (is.matrix(overall_err)) {
    err_means <- colMeans(overall_err, na.rm = TRUE)
    err_sds   <- apply(overall_err, 2, sd, na.rm = TRUE)
    plot(names(err_means), err_means, type = "b", pch = 19,
         ylim = c(0, max(err_means + err_sds, na.rm = TRUE)),
         xlab = "Component", ylab = "Classification Error Rate",
         main = "DIABLO Performance by Component")
    arrows(1:length(err_means), err_means - err_sds,
           1:length(err_means), err_means + err_sds,
           angle = 90, code = 3, length = 0.05)
    abline(h = 1 / length(unique(diablo_res$Y)), lty = 2,
           col = "grey50") # Chance level
  }
  dev.off()

  # 5b. Per-class error rates
  pdf(paste0(out_prefix, "per_class_error.pdf"), width = 8, height = 6)
  if ("class" %in% names(perf_res$error.rate)) {
    for (cl in names(perf_res$error.rate$class)) {
      err_cl <- perf_res$error.rate$class[[cl]]
      if (is.matrix(err_cl)) {
        plot(colMeans(err_cl), type = "b", pch = 19,
             xlab = "Component", ylab = "Error Rate",
             main = paste("Per-Class Error:", cl))
      }
    }
  }
  dev.off()

  # 5c. Selected features summary barplot
  pdf(paste0(out_prefix, "selected_features.pdf"), width = 10, height = 6)
  selected <- selectVar(diablo_res, comp = 1)
  omics_names <- names(selected)

  feat_counts <- sapply(omics_names, function(nm) {
    if (!is.null(selected[[nm]]$value)) nrow(selected[[nm]]$value) else 0
  })

  barplot(feat_counts, col = brewer.pal(length(omics_names), "Set2"),
          main = "Selected Features per Omics — Component 1",
          ylab = "Number of Selected Features",
          las = 2)
  dev.off()

  cat("DIABLO diagnostic plots saved.\n")
}

# ============================================================
# 6. Cross-Omics Feature Pair Network (sPLS/DIABLO)
# ============================================================

plot_cross_omics_network <- function(omics1_features, omics2_features,
                                     cor_threshold = 0.6,
                                     outfile = "outputs/figures/cross_omics_network.pdf") {
  library(igraph)

  # Build nodes
  nodes <- data.frame(
    name = c(omics1_features$feature, omics2_features$feature),
    omics = c(rep("Omics1", nrow(omics1_features)),
              rep("Omics2", nrow(omics2_features))),
    stringsAsFactors = FALSE
  )

  # Build edges (simplified — use actual correlation data)
  # This is a placeholder structure
  edges <- data.frame(
    from = character(0),
    to   = character(0),
    weight = numeric(0)
  )

  if (nrow(edges) == 0) {
    cat("Network: no edges above threshold. Check correlation data.\n")
    return(invisible(NULL))
  }

  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

  # Color nodes by omics type
  V(g)$color <- ifelse(V(g)$omics == "Omics1", "tomato", "steelblue")
  V(g)$size  <- 5

  pdf(outfile, width = 10, height = 10)
  plot(g,
       vertex.label = V(g)$name,
       vertex.label.cex = 0.6,
       vertex.label.dist = 0.5,
       edge.width = abs(E(g)$weight) * 3,
       edge.color = ifelse(E(g)$weight > 0, "red", "blue"),
       main = "Cross-Omics Feature Correlation Network")
  legend("topright",
         legend = c("Omics 1", "Omics 2", "Positive cor", "Negative cor"),
         col = c("tomato", "steelblue", "red", "blue"),
         pch = c(19, 19, NA, NA),
         lty = c(NA, NA, 1, 1),
         cex = 0.8)
  dev.off()

  cat(sprintf("Saved: %s\n", outfile))
}

# ============================================================
# 7. Publication-Ready Multi-Panel Figure
# ============================================================

assemble_multiomics_figure <- function(outfile = "outputs/figures/multiomics_summary.pdf") {
  # This function assumes individual plots have been saved as PDF files
  # and combines them into a multi-panel figure using cowplot

  # Example assembly logic (customize paths as needed):
  # p1 <- ggdraw() + draw_image("outputs/figures/variance_decomposition.pdf")
  # p2 <- ggdraw() + draw_image("outputs/figures/factor_scores.pdf")
  # p3 <- ggdraw() + draw_image("outputs/figures/omics_correlation.pdf")
  #
  # top_row <- plot_grid(p1, p2, ncol = 2, labels = c("A", "B"))
  # final <- plot_grid(top_row, p3, ncol = 1, labels = c("", "C"))
  # ggsave(outfile, final, width = 12, height = 16)

  cat(sprintf("Multi-panel figure: assemble manually using cowplot::plot_grid().\n"))
  cat(sprintf("Output path: %s\n", outfile))
}
```

## Plant-Specific Visualization Notes

- **Subgenome coloring**: For polyploid species (wheat, cotton), color features by subgenome (A/B/D) in heatmaps and network plots. Use distinct palettes: e.g., A=#E41A1C, B=#377EB8, D=#4DAF4A.

- **Chromosome-position plots**: For transcriptome features with genomic coordinates, plot factor weights along chromosomes (similar to Manhattan plot) to identify genomic hotspots of multi-omics coordination.

- **Tissue/Stage faceting**: If samples span tissues or developmental stages, facet sample ordination plots by tissue to reveal tissue-specific factor patterns.

- **Metabolite class annotation**: Use distinct symbols or colors for metabolite classes (flavonoids, terpenoids, alkaloids, lipids, etc.) in loading and network plots.

- **Export formats**: For publication, use `ggsave(..., device = "pdf", width = ..., height = ...)` with vector graphics. For presentations, PNG at 300 DPI.

- **Colorblind-friendly palettes**: Use `RColorBrewer::brewer.pal(..., "Set2")` or `viridis::scale_color_viridis()` throughout.
