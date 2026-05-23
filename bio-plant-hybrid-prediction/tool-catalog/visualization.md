# Hybrid Prediction Visualization

**Goal:** Visualize heterotic groups, GCA/SCA effects, prediction accuracy, and genetic gain
**Best for:** Communicating hybrid breeding results to breeders and stakeholders

## Heterotic Group PCA Plot

```r
library(ggplot2)
library(ggrepel)

# PCA of inbred lines
pca_df <- data.frame(
  PC1 = pca_scores[, 1],
  PC2 = pca_scores[, 2],
  Group = groups$cluster,
  Label = rownames(G)
)

var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(Group))) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  geom_text_repel(aes(label = Label), size = 2.5, max.overlaps = 15) +
  scale_color_brewer(palette = "Set1", name = "Heterotic Group") +
  labs(x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"),
       title = "Heterotic Group Classification by PCA") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("heterotic_groups_pca.png", width = 8, height = 7, dpi = 300)
```

## GCA/SCA Variance Component Plot

```r
# Stacked bar of variance components
var_df <- data.frame(
  Component = c("GCA Male", "GCA Female", "SCA", "Error"),
  Variance = c(gca_male_var, gca_female_var, sca_var, error_var),
  Percentage = c(gca_male_var, gca_female_var, sca_var, error_var) /
               sum(c(gca_male_var, gca_female_var, sca_var, error_var)) * 100
)

ggplot(var_df, aes(x = "", y = Variance, fill = Component)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("#2ecc71", "#3498db", "#e74c3c", "#95a5a6")) +
  labs(title = paste0("Variance Components (Baker's Ratio = ",
                      round(baker_ratio, 3), ")"),
       fill = "Component") +
  theme_void()
```

## GCA Ranking Plot

```r
# Horizontal bar plot of GCA values
gca_rank_df <- data.frame(
  Parent = c(male_names, female_names),
  GCA = c(gca_male_values, gca_female_values),
  Type = c(rep("Male", length(male_names)), rep("Female", length(female_names)))
)
gca_rank_df <- gca_rank_df[order(gca_rank_df$GCA), ]
gca_rank_df$Parent <- factor(gca_rank_df$Parent, levels = gca_rank_df$Parent)

ggplot(gca_rank_df, aes(x = GCA, y = Parent, fill = Type)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Male" = "#3498db", "Female" = "#e74c3c")) +
  labs(x = "GCA (BLUP)", y = "",
       title = "Parent GCA Rankings") +
  theme_minimal()
```

## Prediction Accuracy Scatter Plot

```r
# Observed vs predicted hybrid performance
pred_df <- data.frame(
  Observed = hybrid_data$Yield[test_idx],
  Predicted = predicted_values
)

r_value <- cor(pred_df$Observed, pred_df$Predicted)

ggplot(pred_df, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6, size = 2, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#e74c3c") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = min(pred_df$Observed),
           y = max(pred_df$Predicted) * 0.95,
           label = paste0("r = ", round(r_value, 3)),
           hjust = 0, size = 5) +
  labs(x = "Observed Yield", y = "Predicted Yield",
       title = "Genomic Prediction of Hybrid Performance") +
  theme_minimal()
ggsave("prediction_accuracy.png", width = 7, height = 6, dpi = 300)
```

## Genetic Gain Trajectory

```r
# Projected genetic gain over breeding cycles
cycles <- 1:10
gain <- cumsum(rep(selected_gain_per_cycle, length(cycles)))

ggplot(data.frame(Cycle = cycles, Gain = gain),
       aes(x = Cycle, y = Gain)) +
  geom_line(linewidth = 1.2, color = "#27ae60") +
  geom_point(size = 3, color = "#27ae60") +
  labs(x = "Breeding Cycle", y = "Cumulative Genetic Gain",
       title = "Projected Genetic Gain from Optimal Mate Selection") +
  theme_minimal()
```

## Cross Prediction Heatmap

```r
# Heatmap: predicted performance of all possible crosses
# Rows: Males, Columns: Females
library(pheatmap)

cross_matrix <- matrix(predicted_values,
                       nrow = length(male_ids),
                       ncol = length(female_ids),
                       dimnames = list(male_ids, female_ids))

pheatmap(cross_matrix,
         color = colorRampPalette(c("#3498db", "white", "#e74c3"))(100),
         main = "Predicted Hybrid Performance",
         fontsize_row = 6,
         fontsize_col = 6,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
```

## Plant-Specific Figure Notes

- Multi-panel figure for publication: (a) heterotic group PCA, (b) GCA ranking, (c) prediction accuracy, (d) selected crossing scheme
- Color heterotic groups consistently across all figures (Figure 1 group colors = Figure 2 group colors)
- Include Baker's ratio on variance component plot
- Label top 5-10 parents in GCA ranking plot
- In prediction accuracy plot: show both random CV and forward prediction (if available)

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| PCA points overlap perfectly | Identical lines or low marker density | Increase markers, check for clones/duplicates |
| Heatmap all same color | All predictions similar | Check model convergence |
| Prediction scatter shows clusters | Environment effect not removed | Adjust phenotypes by environment BLUE/BLUP |
| GCA ranking ties visually | Too many parents | Label only top/bottom 10, aggregate middle |
