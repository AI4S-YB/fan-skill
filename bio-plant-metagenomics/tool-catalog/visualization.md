# Metagenomics Visualization

**Goal:** Visualize community composition, functional profiles, MAG quality, and comparative analyses
**Best for:** Publication-ready figures for plant microbiome studies

## Community Composition with phyloseq

```r
library(phyloseq)
library(ggplot2)

# Load data
ps <- readRDS("phyloseq_object.rds")

# Stacked bar plot (phylum level)
ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

plot_bar(ps_rel, fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Relative Abundance") +
  ggtitle("Rhizosphere Community Composition")
```

## PCoA / NMDS Ordination

```r
# Bray-Curtis PCoA
ord <- ordinate(ps, method = "PCoA", distance = "bray")

plot_ordination(ps, ord, color = "Compartment", shape = "Treatment") +
  geom_point(size = 3) +
  stat_ellipse(aes(color = Compartment)) +
  theme_minimal() +
  ggtitle("Beta Diversity by Compartment")

# PERMANOVA test
library(vegan)
dist_matrix <- phyloseq::distance(ps, method = "bray")
adonis2(dist_matrix ~ Compartment * Treatment,
        data = as(sample_data(ps), "data.frame"))
```

## MAG Quality Visualization

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load CheckM output
qc = pd.read_csv("mag_quality_report.tsv", sep="\t")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Histogram of completeness
axes[0].hist(qc["Completeness"], bins=30, color="#3498db", edgecolor="black")
axes[0].axvline(x=90, color="green", linestyle="--", label="High quality (90%)")
axes[0].axvline(x=50, color="orange", linestyle="--", label="Medium quality (50%)")
axes[0].set_xlabel("Completeness (%)")
axes[0].set_ylabel("Number of MAGs")
axes[0].legend()
axes[0].set_title("MAG Completeness Distribution")

# Scatter: completeness vs contamination
axes[1].scatter(qc["Completeness"], qc["Contamination"],
                alpha=0.6, c="#e74c3c", edgecolors="black", linewidth=0.5)
axes[1].axhline(y=5, color="green", linestyle="--", label="Contamination < 5%")
axes[1].axvline(x=90, color="green", linestyle="--")
axes[1].axvline(x=50, color="orange", linestyle="--")
axes[1].set_xlabel("Completeness (%)")
axes[1].set_ylabel("Contamination (%)")
axes[1].legend()
axes[1].set_title("MAG Quality Overview")

plt.tight_layout()
plt.savefig("mag_quality.png", dpi=300)
```

## Functional Heatmap

```python
import seaborn as sns

# KEGG module completeness matrix
modules = pd.read_csv("kegg_modules.csv", index_col=0)

plt.figure(figsize=(14, 10))
sns.clustermap(modules,
               cmap="RdYlBu_r",
               method="ward",
               figsize=(14, 12),
               xticklabels=True,
               cbar_kws={'label': 'Module Completeness (%)'})
plt.title("KEGG Module Completeness Across MAGs")
plt.savefig("kegg_module_heatmap.png", dpi=300, bbox_inches='tight')
```

## Alpha Diversity by Compartment

```r
# Shannon diversity
alpha <- estimate_richness(ps, measures = "Shannon")
alpha$Sample <- rownames(alpha)
alpha <- merge(alpha, sample_data(ps), by = "row.names")

ggplot(alpha, aes(x = Compartment, y = Shannon, fill = Compartment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_minimal() +
  ylab("Shannon Diversity Index") +
  ggtitle("Alpha Diversity by Plant Compartment")
```

## Plant-Specific Figure Notes

- Multi-panel figure: (a) community bar plot, (b) PCoA, (c) alpha diversity, (d) MAG quality scatter
- Color by plant compartment (rhizosphere red, endosphere blue, phyllosphere green, soil brown)
- Always show relative abundance (%), not raw counts
- For endosphere: indicate plant DNA removal efficiency in supplementary figure

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| PCoA shows no separation | Similar communities or too few samples | Add PERMANOVA statistics |
| Heatmap too dense | Too many MAGs | Filter to quality MAGs only or aggregate by genus |
| Phyloseq import fails | Incompatible OTU/tax table formats | Use `import_biom()` or check column names |
| Plot colors indistinguishable | Too many categories | Aggregate rare taxa to "Other" |
