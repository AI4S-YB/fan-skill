# Diversity Analysis

**Goal:** Calculate alpha and beta diversity metrics and test for group differences
**Best for:** All amplicon studies — core analysis for comparing microbial communities

## Prerequisites
- QIIME 2 with diversity plugin
- Feature table (QZA), rooted phylogenetic tree (QZA), metadata (TSV)
- Rarefied table (recommended for alpha diversity)

## Alpha Diversity

### Calculate Alpha Diversity Metrics

```bash
qiime diversity alpha \
  --i-table rarefied-table.qza \
  --p-metric shannon \
  --o-alpha-diversity shannon.qza

qiime diversity alpha \
  --i-table rarefied-table.qza \
  --p-metric observed_features \
  --o-alpha-diversity observed_features.qza

qiime diversity alpha \
  --i-table rarefied-table.qza \
  --p-metric faith_pd \
  --i-phylogeny rooted-tree.qza \
  --o-alpha-diversity faith_pd.qza

qiime diversity alpha \
  --i-table rarefied-table.qza \
  --p-metric pielou_e \
  --o-alpha-diversity evenness.qza
```

### Alpha Diversity Group Comparison

```bash
# Kruskal-Wallis test (non-parametric)
qiime diversity alpha-group-significance \
  --i-alpha-diversity shannon.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization shannon-group-significance.qzv

# Alpha diversity boxplot
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 10000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```

## Beta Diversity

### Calculate Beta Diversity Matrices

```bash
# Bray-Curtis (abundance-weighted, non-phylogenetic)
qiime diversity beta \
  --i-table rarefied-table.qza \
  --p-metric braycurtis \
  --o-distance-matrix braycurtis.qza

# Weighted UniFrac (abundance-weighted, phylogenetic)
qiime diversity beta-phylogenetic \
  --i-table rarefied-table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-metric weighted_unifrac \
  --o-distance-matrix weighted_unifrac.qza

# Unweighted UniFrac (presence/absence, phylogenetic)
qiime diversity beta-phylogenetic \
  --i-table rarefied-table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-metric unweighted_unifrac \
  --o-distance-matrix unweighted_unifrac.qza

# Jaccard (presence/absence, non-phylogenetic)
qiime diversity beta \
  --i-table rarefied-table.qza \
  --p-metric jaccard \
  --o-distance-matrix jaccard.qza
```

### PCoA Ordination

```bash
qiime diversity pcoa \
  --i-distance-matrix braycurtis.qza \
  --o-pcoa braycurtis-pcoa.qza

qiime diversity pcoa \
  --i-distance-matrix weighted_unifrac.qza \
  --o-pcoa weighted-unifrac-pcoa.qza
```

### Beta Diversity Group Significance (PERMANOVA)

```bash
qiime diversity beta-group-significance \
  --i-distance-matrix braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Treatment \
  --p-method permanova \
  --p-pairwise \
  --o-visualization braycurtis-permanova.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix weighted_unifrac.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Tissue \
  --p-method permanova \
  --p-pairwise \
  --o-visualization unifrac-permanova.qzv
```

### PERMDISP (Dispersion Test)

```bash
qiime diversity beta-group-significance \
  --i-distance-matrix braycurtis.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Treatment \
  --p-method permdisp \
  --o-visualization dispersion-test.qzv
```

## Key Metrics Overview

| Metric | Type | Weighting | Phylogenetic | Use Case |
|--------|------|-----------|--------------|----------|
| Shannon | Alpha | Abundance | No | General diversity |
| Observed features | Alpha | Presence/Absence | No | Richness |
| Faith's PD | Alpha | Presence/Absence | Yes | Phylogenetic richness |
| Pielou's evenness | Alpha | Abundance | No | Evenness |
| Bray-Curtis | Beta | Abundance | No | Community composition |
| Weighted UniFrac | Beta | Abundance | Yes | Phylogenetic composition |
| Unweighted UniFrac | Beta | Presence/Absence | Yes | Phylogenetic membership |
| Jaccard | Beta | Presence/Absence | No | Community membership |

## Rarefaction Depth Selection

Check the feature table summary to select sampling depth:
- Choose depth that retains > 90% of samples
- Typical depths: 5000-10000 (soil), 2000-5000 (rhizosphere), 1000-3000 (endosphere)

## Plant-Specific Considerations

- Rhizosphere vs bulk soil comparisons should use both Bray-Curtis and UniFrac to see phylogenetic signal
- For plant genotype comparisons, Faith's PD is informative — it captures phylogenetic breadth of the recruited microbiome
- Root compartment (rhizosphere vs rhizoplane vs endosphere) is typically the strongest factor in beta diversity; include it as a fixed effect
- Plant developmental stage affects alpha diversity — include it as a covariate when samples span multiple time points

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Distance matrix contains NaN" | Sample has zero reads after rarefaction | Drop empty samples or lower rarefaction depth |
| "Phylogenetic tree required" | Using phylogenetic metric without tree | Build tree with `qiime phylogeny align-to-tree-mafft-fasttree` |
| PERMANOVA not significant | Low effect size or small sample size | Check dispersion (PERMDISP); consider increasing replicates |
| Faith's PD values suspicious | Tree poorly resolved | Check tree quality; use more reference sequences |
