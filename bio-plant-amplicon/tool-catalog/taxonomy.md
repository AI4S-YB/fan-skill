# Taxonomy Assignment

**Goal:** Assign taxonomic classifications to ASV/OTU representative sequences
**Best for:** All amplicon studies — critical for biological interpretation

## Prerequisites
- QIIME 2 with feature-classifier plugin
- Pre-trained classifier or reference database
- Representative sequences (FASTA or QZA)

## SILVA 138 for 16S rRNA

### Download and Train Classifier

```bash
# Download SILVA 138 99% OTUs (16S only)
wget https://data.qiime2.org/2024.2/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2024.2/common/silva-138-99-tax.qza

# Extract V4 region for classifier training (if using V4 primers)
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 250 \
  --o-reads silva-138-99-v4-ref-seqs.qza

# Train classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-v4-ref-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-v4-classifier.qza
```

### Assign Taxonomy

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-v4-classifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence 0.7 \
  --o-classification taxonomy.qza

# Visualize taxonomy
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Generate taxonomy bar plot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-barplot.qzv
```

## UNITE for ITS (Fungi)

```bash
# Download UNITE QIIME 2 classifier (pre-trained)
wget https://files.plutof.ut.ee/qiime2/classifiers/unite-ver9-classifier-2024.qza

# Assign taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier unite-ver9-classifier-2024.qza \
  --i-reads rep-seqs.qza \
  --p-confidence 0.7 \
  --o-classification taxonomy.qza
```

## Alternative: BLAST-based Classification

```bash
# For custom databases or non-standard markers
qiime feature-classifier classify-consensus-blast \
  --i-query rep-seqs.qza \
  --i-reference-reads custom-ref-seqs.qza \
  --i-reference-taxonomy custom-ref-tax.qza \
  --p-maxaccepts 1 \
  --p-perc-identity 0.97 \
  --o-classification blast-taxonomy.qza
```

## Filter Plant Contaminants

```bash
# Filter out chloroplast and mitochondrial reads
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude "Chloroplast,Mitochondria" \
  --o-filtered-table table-no-plant.qza

# Also remove unassigned sequences
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include p__ \
  --o-filtered-table table-classified.qza
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| --p-confidence | Classification confidence threshold (0.7 for genus, 0.8 for phylum) |
| --p-f-primer | Forward primer sequence for read extraction |
| --p-r-primer | Reverse primer sequence |
| --p-trunc-len | Expected amplicon length after trimming |
| --p-maxaccepts | Max hits to consider (BLAST) |
| --p-perc-identity | Minimum percent identity (BLAST) |
| --p-exclude | Taxa to exclude (comma-separated) |

## Plant-Specific Considerations

- ALWAYS filter Chloroplast and Mitochondria after taxonomy assignment for plant samples
- For nifH, amoA, or other functional gene amplicons, use custom databases (e.g., FunGene, custom BLAST DBs)
- Root endosphere samples: expect 30-60% chloroplast/mitochondria reads even with PNA clamps
- Soil/rhizosphere: expect < 5% plant-derived reads

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Classifier training failed" | Insufficient memory | Use pre-trained classifier or reduce reference size |
| Most ASVs "Unassigned" | Wrong primer region or database | Verify primer sequences match the extracted region |
| "No reads after taxonomy filter" | All sequences classified as chloroplast | Check if PNA clamp was used; verify primer specificity |
| Confidence too low for genus | Short amplicon or conserved region | Lower confidence threshold for genus-level calls |
