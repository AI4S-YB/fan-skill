# QIIME 2 Workflow

**Goal:** End-to-end amplicon analysis pipeline — import, denoise, classify, analyze diversity
**Best for:** Standardized amplicon analysis with reproducible provenance tracking

## Prerequisites
- QIIME 2 2024.2+
- Conda environment with qiime2
- Metadata file (TSV format) with sample IDs and group information

## Import Data

```bash
# Import paired-end FASTQ (Casava 1.8 format)
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path raw_reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux.qza

# Visualize quality
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux-summary.qzv
```

## Denoising Options

### DADA2 (Recommended)

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 240 --p-trunc-len-r 200 \
  --p-max-ee-f 2 --p-max-ee-r 2 \
  --p-n-threads 16 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza
```

### Deblur (Fixed-Length ASV)

```bash
# Deblur requires single-end reads; first join paired reads
qiime vsearch join-pairs \
  --i-demultiplexed-seqs demux.qza \
  --o-joined-seqs joined.qza

# Quality filter
qiime quality-filter q-score \
  --i-demux joined.qza \
  --o-filtered-sequences filtered.qza \
  --o-filter-stats filter-stats.qza

# Deblur denoising
qiime deblur denoise-16S \
  --i-demultiplexed-seqs filtered.qza \
  --p-trim-length 250 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-stats deblur-stats.qza
```

## Feature Table Operations

```bash
# Summarize feature table
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table-summary.qzv \
  --m-sample-metadata-file metadata.tsv

# Filter features by frequency (remove rare ASVs)
qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-frequency 10 \
  --p-min-samples 2 \
  --o-filtered-table filtered-table.qza

# Rarefy to even sampling depth
qiime feature-table rarefy \
  --i-table filtered-table.qza \
  --p-sampling-depth 5000 \
  --o-rarefied-table rarefied-table.qza
```

## Phylogenetic Tree Construction

```bash
# Align sequences
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| --p-trunc-len-f/r | Trim forward/reverse reads to this length |
| --p-max-ee | Maximum expected errors per read |
| --p-trim-length | Trim all reads to fixed length (Deblur) |
| --p-min-frequency | Minimum feature frequency filter |
| --p-min-samples | Minimum sample occurrence filter |
| --p-sampling-depth | Rarefaction depth |

## Plant-Specific Considerations

- Use `--p-min-frequency` to remove likely contaminant ASVs in low-biomass samples
- For rhizosphere samples, 5000-10000 reads/sample is common
- For endosphere/phyllosphere, 1000-3000 reads/sample may be necessary; note low depth in methods
- Always include extraction blank and PCR negative controls; use their profiles to filter contaminants

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Plugin error from dada2" | Installation issue or wrong QIIME version | Reinstall qiime2 in clean conda env |
| "Metadata column not found" | Column name mismatch | Verify with `head metadata.tsv` |
| "No sequences after rarefaction" | Sampling depth too high | Check `table-summary.qzv` for per-sample counts |
| "Casava format error" | Wrong FASTQ naming | Use `CasavaOneEightSingleLanePerSampleDirFmt` or rename files |
