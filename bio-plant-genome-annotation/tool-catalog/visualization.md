# Annotation Visualization

**Goal:** Generate publication-quality summary plots and genome browser tracks for plant genome annotations
**Best for:** Visualizing annotation results for publications and presentations

## Prerequisites
- R 4.0+
- R packages: ggplot2, circlize, GenomicRanges, rtracklayer, ComplexHeatmap
- Python 3 with matplotlib, seaborn (alternative)

## Annotation Summary Dashboard

```r
library(ggplot2)
library(reshape2)

# Read annotation statistics
stats <- data.frame(
  Category = c("Genes", "mRNA", "CDS", "Exons", "Introns"),
  Count = c(45231, 67219, 321987, 256768, 254769)
)

# Gene distribution across chromosomes
chr_stats <- read.table("gene_count_per_chr.txt", header = TRUE)

pdf("outputs/figures/annotation_summary.pdf", width = 14, height = 10)

# Gene count per chromosome
ggplot(chr_stats, aes(x = reorder(Chr, -GeneCount), y = GeneCount)) +
  geom_bar(stat = "identity", fill = "#2E86AB") +
  labs(title = "Gene Distribution Across Chromosomes",
       x = "Chromosome", y = "Number of Genes") +
  theme_minimal(base_size = 14)
```

## BUSCO Visualization

```r
# Create BUSCO summary plot from short_summary.txt
busco_data <- data.frame(
  Category = c("Complete (S)", "Complete (D)", "Fragmented", "Missing"),
  Count = c(1407, 53, 82, 72),
  Color = c("#4CB944", "#2E86AB", "#F18F01", "#C73E1D")
)

ggplot(busco_data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = busco_data$Color) +
  labs(title = "BUSCO Assessment (embryophyta_odb10, n=1614)") +
  theme_void()
```

## Genome Browser Track Setup

### Convert GFF3 to BigBed/BigWig

```bash
# Sort GFF3
gff3sort.pl predicted_genes.gff3 > predicted_genes_sorted.gff3

# Convert to BED
gff2bed < predicted_genes_sorted.gff3 > predicted_genes.bed

# Convert to BigBed (for UCSC/JBrowse)
bedToBigBed predicted_genes.bed chrom.sizes predicted_genes.bb

# Generate coverage track (if you have alignments)
bamCoverage -b rnaseq_sorted.bam -o rnaseq_coverage.bw -p 8 \
  --normalizeUsing CPM --binSize 10
```

### JBrowse2 Configuration

```json
{
  "tracks": [
    {
      "type": "FeatureTrack",
      "trackId": "genes",
      "name": "Predicted Genes",
      "assemblyNames": ["species_assembly"],
      "adapter": {
        "type": "Gff3TabixAdapter",
        "gffGzLocation": { "uri": "predicted_genes_sorted.gff3.gz" }
      }
    },
    {
      "type": "QuantitativeTrack",
      "trackId": "rnaseq_cov",
      "name": "RNA-seq Coverage (CPM)",
      "adapter": {
        "type": "BigWigAdapter",
        "bigWigLocation": { "uri": "rnaseq_coverage.bw" }
      }
    }
  ]
}
```

## Repeat Landscape Plot

```r
library(ggplot2)

# Read RepeatMasker output
repeats <- read.table("genome.fasta.out", skip = 3, fill = TRUE)

# Kimura distance distribution (repeat age estimation)
pdf("outputs/figures/repeat_landscape.pdf", width = 10, height = 6)
ggplot(repeats, aes(x = V2, fill = V11)) +
  geom_histogram(binwidth = 1, position = "stack") +
  labs(title = "Repeat Landscape (Kimura Divergence)",
       x = "Kimura substitution level (%)",
       y = "Genome proportion (%)") +
  scale_fill_discrete(name = "Repeat Class") +
  theme_minimal(base_size = 12)
dev.off()
```

## Gene Feature Length Distribution

```r
# Plot distributions: gene length, exon count, CDS length, intron length
features <- read.table("gene_features.txt", header = TRUE)

# Multi-panel distribution plot
pdf("outputs/figures/gene_feature_distribution.pdf", width = 12, height = 10)

p1 <- ggplot(features, aes(x = gene_length)) +
  geom_density(fill = "#2E86AB", alpha = 0.6) +
  labs(x = "Gene Length (bp)", y = "Density") + theme_minimal()

p2 <- ggplot(features, aes(x = exon_count)) +
  geom_bar(fill = "#4CB944", alpha = 0.7) +
  labs(x = "Exon Count per Gene", y = "Frequency") + theme_minimal()

p3 <- ggplot(features, aes(x = cds_length)) +
  geom_density(fill = "#F18F01", alpha = 0.6) +
  labs(x = "CDS Length (bp)", y = "Density") + theme_minimal()

p4 <- ggplot(features, aes(x = intron_length)) +
  geom_density(fill = "#C73E1D", alpha = 0.6) +
  scale_x_log10() +
  labs(x = "Intron Length (bp, log10)", y = "Density") + theme_minimal()

grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| --normalizeUsing | Coverage normalization method (CPM, RPKM, BPM) |
| --binSize | Resolution for coverage tracks |
| binwidth | Histogram bin width for Kimura distance |
| chrom.sizes | Chromosome size file for BigBed conversion |

## Plant-Specific Considerations

- For polyploid genomes, color-code genes by subgenome in tracks
- Include centromere/pericentromere annotation if available
- Repeat landscape plots are particularly informative for plant genomes (LTR burst history)
- BUSCO pie charts are the standard for plant genome paper figures

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "gff2bed not found" | BEDOPS not installed | Install BEDOPS or use `gffread` conversion |
| chrom.sizes mismatch | Wrong assembly version | Regenerate with `samtools faidx` |
| Empty tracks in JBrowse | Compression or indexing missing | Ensure .tbi index for VCF/GFF, .bai for BAM |
| R memory exhausted | Large genome GFF loaded entirely | Use GenomicFiles for streaming or sample chromosomes |
