# Microbiome Visualization

**Goal:** Generate publication-quality figures for plant microbiome amplicon studies
**Best for:** All amplicon studies — taxonomy bar plots, PCoA ordination, heatmaps, and network plots

## Prerequisites
- R 4.0+
- R packages: ggplot2, phyloseq, vegan, pheatmap, RColorBrewer, igraph, ggraph
- QIIME 2 (for built-in visualization)

## Taxonomy Bar Plot

### QIIME 2 Version

```bash
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-barplot.qzv
```

### Custom R Version (Better Control)

```r
library(phyloseq)
library(ggplot2)
library(dplyr)

# Aggregate to phylum level
physeq_phylum <- tax_glom(physeq, taxrank = "Phylum")
# Transform to relative abundance
physeq_rel <- transform_sample_counts(physeq_phylum, function(x) x / sum(x) * 100)

# Melt for ggplot2
melted <- psmelt(physeq_rel)

# Group low-abundance phyla as "Other"
top_phyla <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE))[1:8]
melted$Phylum[!melted$OTU %in% top_phyla] <- "Other"

# Plot
ggplot(melted, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Treatment, scales = "free_x", space = "free") +
  labs(y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set3")
```

## PCoA Ordination Plot

### QIIME 2 Version

```bash
# Emperor plot (interactive 3D)
qiime emperor plot \
  --i-pcoa braycurtis-pcoa.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization braycurtis-emperor.qzv
```

### Custom R Version

```r
library(phyloseq)
library(ggplot2)
library(vegan)

# Calculate Bray-Curtis distance and PCoA
dist_bc <- distance(physeq_rel, method = "bray")
pcoa <- ordinate(physeq_rel, method = "PCoA", distance = dist_bc)

# Extract variance explained
eigenvals <- pcoa$values$Relative_eig
var1 <- round(eigenvals[1] * 100, 1)
var2 <- round(eigenvals[2] * 100, 1)

# Plot
plot_ordination(physeq_rel, pcoa, color = "Treatment", shape = "Tissue") +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = Treatment), level = 0.95) +
  labs(title = "Bray-Curtis PCoA",
       x = paste0("PCoA 1 (", var1, "%)"),
       y = paste0("PCoA 2 (", var2, "%)")) +
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Set1")
```

## Alpha Diversity Boxplot

```r
# Extract alpha diversity metrics
alpha_div <- estimate_richness(physeq, measures = c("Shannon", "Observed", "Faith"))
alpha_div$Sample <- rownames(alpha_div)
alpha_div <- merge(alpha_div, sample_data_df, by.x = "Sample", by.y = "row.names")

# Boxplot
ggplot(alpha_div, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  stat_compare_means(method = "kruskal.test", label.y = max(alpha_div$Shannon) * 1.1) +
  labs(y = "Shannon Diversity Index", x = "") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("#2E86AB", "#F18F01", "#C73E1D"))
```

## Heatmap of Differential Taxa

```r
library(pheatmap)

# Extract significant genera from ANCOM-BC
sig_genera <- sig_taxa$taxon
physeq_sig <- prune_taxa(sig_genera, physeq_genus)

# CLR transform
physeq_clr <- microbiome::transform(physeq_sig, "clr")

# Extract matrix
mat <- as.matrix(otu_table(physeq_clr))

# Row annotation by treatment
annotation_col <- data.frame(
  Treatment = sample_data(physeq_sig)$Treatment,
  Tissue = sample_data(physeq_sig)$Tissue
)
rownames(annotation_col) <- sample_names(physeq_sig)

pheatmap(mat,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("#2E86AB", "white", "#C73E1D"))(100),
         fontsize_row = 8,
         main = "Differentially Abundant Genera (CLR)")
```

## Co-occurrence Network

```r
library(igraph)
library(ggraph)

# Calculate correlation matrix
counts <- as.matrix(otu_table(physeq_genus))
cor_mat <- cor(t(counts), method = "spearman")

# Keep significant correlations only (|rho| > 0.6, p < 0.05)
cor_mat[abs(cor_mat) < 0.6] <- 0
diag(cor_mat) <- 0

# Build graph
graph <- graph_from_adjacency_matrix(cor_mat, mode = "undirected", weighted = TRUE)

# Plot
ggraph(graph, layout = "fr") +
  geom_edge_link(aes(edge_width = abs(weight), edge_alpha = abs(weight)),
                 color = "grey50") +
  geom_node_point(aes(size = degree(graph)), color = "#2E86AB") +
  theme_void() +
  labs(title = "Genus Co-occurrence Network (|rho| > 0.6)")
```

## Plant-Specific Visualization Notes

- Compartment effect (rhizosphere vs endosphere) is often the dominant PCoA axis — use shape or color to show it
- For plant genotype studies, use faceted bar plots to show genotype-specific enrichment
- Soil microbiome studies with spatial sampling: add spatial coordinates to PCoA to check for distance-decay patterns
- Use plant-friendly color palettes (green-to-brown for compartments, or crop-specific colors)

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| stat_ellipse level | Confidence level for ellipse (0.95 = 95% CI) |
| scale_fill_brewer palette | ColorBrewer palette name |
| clustering_distance_rows | Distance metric for heatmap clustering |
| cor method | Correlation method (spearman for non-normal) |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| PCoA plot crowded | Too many samples | Use faceting or transparency |
| Heatmap all same color | Not scaled properly | Use `scale = "row"` for relative differences |
| Emperor blank | Missing metadata | Verify metadata column names match |
| Co-occurrence graph too dense | Correlation threshold too low | Increase |rho| threshold |
