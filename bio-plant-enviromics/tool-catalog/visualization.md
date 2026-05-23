# Enviromics Visualization (环境组学可视化)

**Goal:** Visualize environmental characterization, G×E patterns, stability,
and TPE analyses

## Prerequisites
- R 4.0+, packages: `ggplot2`, `metan`, `factoextra`, `corrplot`, `pheatmap`,
  `patchwork`, `sf`, `maps`
- Outputs from env-index, gxe-modeling, tpe-analysis, stability

## Visualization Catalog

### 1. Environmental Index Visualization

```r
library(ggplot2)
library(factoextra)
library(patchwork)

# PCA scree plot
p1 <- fviz_eig(pca_res, addlabels = TRUE,
               main = "Variance Explained by Environmental PCs")

# PCA variable contributions
p2 <- fviz_pca_var(pca_res,
                    col.var = "contrib",
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE,
                    title = "Climate Variable Contributions")

# PCA biplot (environments + variables)
p3 <- fviz_pca_biplot(pca_res,
                       repel = TRUE,
                       col.var = "#2E9FDF",
                       col.ind = "#E7B800",
                       title = "Environmental PCA Biplot")

# Environmental index bar plot
ggplot(env_index, aes(x = reorder(Environment, PC1), y = PC1, fill = PC1)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  labs(x = "Environment", y = "Environmental Index (PC1)",
       title = "Environmental Gradient") +
  coord_flip() + theme_minimal()
```

### 2. GGE/AMMI Biplots

```r
library(metan)

# GGE biplot - comprehensive visualization
gge_model <- gge(data, env = Environment, gen = Genotype, resp = Yield)

# Which-Won-Where pattern
plot(gge_model, type = 3,size.shape.winner.line = 1,
     col.gen = "blue", col.env = "red",
     title = "Which-Won-Where: GGE Biplot")

# Mean vs Stability
plot(gge_model, type = 2,
     col.gen = "blue", col.env = "red",
     title = "Mean Performance vs Stability")

# Ranking genotypes relative to ideal
plot(gge_model, type = 8,
     title = "Genotype Ranking (Ideal Environment)")

# Environment relationship
plot(gge_model, type = 6,
     title = "Environment Relationship Biplot")
```

### 3. Genetic Correlation Heatmap (from FA model)

```r
library(pheatmap)
library(corrplot)

# Genetic correlation matrix between environments
# gen_cor_mat: from FA model output

pheatmap(gen_cor_mat,
         color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "Genetic Correlations Between Environments",
         fontsize_number = 8,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# Correlogram alternative
corrplot(gen_cor_mat,
         method = "color",
         type = "upper",
         addCoef.col = "black",
         number.cex = 0.7,
         tl.col = "black",
         tl.cex = 0.8,
         title = "Genetic Correlations",
         mar = c(0, 0, 2, 0))
```

### 4. Reaction Norm Plot

```r
# Reaction norms: genotype response to environmental gradient
ggplot(dat, aes(x = PC1, y = Yield, group = Genotype, color = Genotype)) +
  stat_smooth(method = "lm", se = FALSE, alpha = 0.3, linewidth = 0.5) +
  geom_point(data = dat %>%
               group_by(Genotype, Environment) %>%
               summarise(Yield = mean(Yield), PC1 = mean(PC1), .groups = "drop"),
             alpha = 0.6, size = 2) +
  labs(x = "Environmental Index (PC1)",
       y = "Grain Yield",
       title = "Reaction Norms: Genotype Response to Environment") +
  theme_minimal() +
  theme(legend.position = "none")
```

### 5. Stability Plots

```r
# Mean vs CV (static stability)
stability_res$stats %>%
  ggplot(aes(x = Y, y = CV, label = GEN)) +
  geom_point(aes(color = CV), size = 2) +
  geom_text(vjust = -0.8, size = 3) +
  geom_vline(xintercept = mean(stability_res$stats$Y), linetype = "dashed", alpha = 0.5) +
  scale_color_gradient(low = "green", high = "red") +
  labs(x = "Mean Yield", y = "Coefficient of Variation (%)",
       title = "Static Stability: Mean vs CV") +
  theme_minimal()

# FW regression: slope vs mean yield
stability_res$stats %>%
  ggplot(aes(x = Y, y = FW, label = GEN)) +
  geom_point(aes(color = abs(FW - 1)), size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  geom_text(vjust = -0.8, size = 3) +
  labs(x = "Mean Yield", y = "Finlay-Wilkinson Slope (b)",
       title = "Dynamic Stability: Regression Slope vs Mean") +
  theme_minimal()
```

### 6. TPE Visualization

```r
# TPE clusters on PCA plot
env_data$TPE <- paste0("TPE_", km$cluster)

ggplot(env_data, aes(x = PC1, y = PC2, color = TPE, label = rownames(env_data))) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.68) +
  geom_text(vjust = -0.8, size = 3) +
  labs(x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
       title = "Target Population of Environments (TPE)") +
  theme_minimal()

# TPE dendrogram
fviz_dend(hc, k = optimal_k,
          rect = TRUE, rect_fill = TRUE,
          main = "TPE Hierarchical Clustering")
```

### 7. Map Visualization of Trial Locations

```r
library(sf)
library(maps)
library(rnaturalearth)

# Trial location map with environmental overlay
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70") +
  geom_point(data = locations,
             aes(x = lon, y = lat, color = PC1_score, size = Yield),
             alpha = 0.8) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue") +
  coord_sf(xlim = range(locations$lon) + c(-5, 5),
           ylim = range(locations$lat) + c(-5, 5)) +
  labs(title = "Trial Locations with Environmental Index",
       color = "Env. Index (PC1)", size = "Mean Yield") +
  theme_minimal()
```

## Plant Relevance

- **Which-Won-Where is the most actionable plot for breeders**: It directly
  identifies which genotype(s) perform best in which mega-environment.
- **Color choices matter**: Use red-blue gradients for temperature (red=hot,
  blue=cool), green-brown for moisture (green=wet, brown=dry). Don't use
  red-green scales (colorblind-unfriendly).
- **Stability biplots**: Always label the best and worst genotypes. Breeders
  want to know names, not just points.
- **Publication quality**: For journal figures, use `theme_classic()` or
  `theme_bw()`, remove unnecessary gridlines, and use at least 10pt font.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| Overlapping labels in biplot | Too many genotypes/environments | Use `ggrepel` or show only top/bottom performers |
| Correlation heatmap all 1.0 | Using raw data instead of correlation matrix | Verify FA model output before plotting |
| Map projection distortion | Wrong CRS for region | Use appropriate UTM zone or regional projection |
