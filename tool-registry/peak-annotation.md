# Peak Annotation: ChIP-seq/ATAC-seq Peak Functional Annotation

## Tool Overview

**Tool ID**: `peak-annotation`
**Category**: Epigenomics Analysis
**Purpose**: Annotate genomic peaks to genes and regulatory elements
**Primary Tool**: ChIPseeker (R/Bioconductor)
**Alternative**: HOMER annotatePeaks.pl, BEDTools

## When to Use

- ChIP-seq peak annotation to nearest genes
- ATAC-seq peak functional annotation
- Promoter/enhancer identification
- Peak enrichment analysis by genomic features

## Key Parameters

| Parameter | Default | Recommended for Plants | Description |
|-----------|---------|------------------------|-------------|
| `tssRegion` | c(-3000, 3000) | **c(-2000, 500)** | TSS region for promoter definition |
| `level` | "gene" | "gene" | Annotation level (gene/transcript) |
| `assignGenomicAnnotation` | TRUE | TRUE | Assign genomic annotation |
| `annoDb` | NULL | Species-specific | Organism annotation database |

### Critical Parameter: `tssRegion`

**Why c(-2000, 500) for plants?**

| Species Type | Recommended tssRegion | Rationale |
|--------------|----------------------|-----------|
| **Plants** | c(-2000, 500) | Plant core promoters are compact; TATA box typically -30 to -100 |
| Mammals | c(-3000, 3000) | Mammalian promoters can be more extended |
| Yeast | c(-500, 200) | Very compact promoters |

**Plant-specific reasoning**:
- Plant promoters typically have smaller core regions
- Upstream regulatory elements often within 2kb
- 5' UTR shorter in plants (median ~150bp vs ~300bp in mammals)

## Recommended Parameter Sets

### Plant ChIP-seq Annotation

```r
library(ChIPseeker)
library(GenomicFeatures)

# Load plant annotation (create from GFF3 if not available)
txdb <- makeTxDbFromGFF("species_annotation.gff3", format="GFF3")

# Annotate peaks
peak_anno <- annotatePeak(
  peak = "peaks.narrowPeak",
  TxDb = txdb,
  tssRegion = c(-2000, 500),
  level = "gene",
  assignGenomicAnnotation = TRUE,
  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")
)
```

### ATAC-seq Peak Annotation

```r
# For ATAC-seq, use broader TSS region
peak_anno <- annotatePeak(
  peak = "atac_peaks.narrowPeak",
  TxDb = txdb,
  tssRegion = c(-2000, 1000),  # ATAC peaks can be in broader regulatory regions
  level = "gene"
)
```

## Annotation Categories

| Category | Description | Typical % in Plants |
|----------|-------------|---------------------|
| **Promoter** | Within tssRegion upstream of TSS | 10-30% |
| **5' UTR** | 5' untranslated region | 2-5% |
| **3' UTR** | 3' untranslated region | 2-5% |
| **Exon** | Coding exon | 5-15% |
| **Intron** | Intronic region | 20-40% |
| **Downstream** | Downstream of gene end | 1-5% |
| **Intergenic** | Between genes | 20-50% |

## Output Files

| Output | Content |
|--------|---------|
| Annotation table | Peak ID, genomic location, annotation, gene ID, distance to TSS |
| Pie chart | Distribution of peaks across genomic features |
| Upset plot | Overlap of peaks with multiple features |
| Heatmap | Peak signal around TSS |

## Plant-Specific Considerations

### 1. Creating TxDb from GFF3

Most plant species don't have pre-built TxDb packages:

```r
# Create TxDb from GFF3
txdb <- makeTxDbFromGFF(
  file = "species_annotation.gff3",
  format = "GFF3",
  organism = "Species name"
)

# Save for future use
saveDb(txdb, file="species_txdb.sqlite")
```

### 2. Handling Polyploid Genomes

For polyploid species with subgenomes:

```r
# Annotate peaks separately for each subgenome
# Split peaks by subgenome assignment
subgenome_a_peaks <- peaks[seqnames(peaks) %in% chrA]
subgenome_b_peaks <- peaks[seqnames(peaks) %in% chrB]

# Annotate each set
anno_a <- annotatePeak(subgenome_a_peaks, TxDb = txdb, tssRegion = c(-2000, 500))
anno_b <- annotatePeak(subgenome_b_peaks, TxDb = txdb, tssRegion = c(-2000, 500))

# Compare annotation distributions
```

### 3. Gene Density Bias

Plant genomes have variable gene density:

```r
# Account for gene density in enrichment analysis
# Use random peaks as background matched for gene density
background <- createRandomPeaks(
  peaks,
  genome = genome,
  match.geneDensity = TRUE
)
```

## Integration with ChIP-seq Pipeline

```r
# Complete annotation workflow

library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)

# 1. Load peaks
peaks <- readPeakFile("peaks.narrowPeak")

# 2. Create/load TxDb
txdb <- loadDb("species_txdb.sqlite")

# 3. Annotate peaks
peak_anno <- annotatePeak(
  peak = peaks,
  TxDb = txdb,
  tssRegion = c(-2000, 500),
  level = "gene"
)

# 4. Get annotation data frame
anno_df <- as.data.frame(peak_anno)

# 5. Functional enrichment of target genes
target_genes <- unique(anno_df$geneId)
ego <- enrichGO(
  gene = target_genes,
  OrgDb = orgdb,  # Use appropriate organism database
  ont = "BP",
  pAdjustMethod = "BH"
)

# 6. Visualize
plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno)
```

## Functional Enrichment Analysis

After annotation, perform GO/KEGG enrichment:

```r
# Get genes with peaks in promoters
promoter_genes <- anno_df %>%
  filter(annotation == "Promoter") %>%
  pull(geneId) %>%
  unique()

# GO enrichment
ego <- enrichGO(
  gene = promoter_genes,
  OrgDb = org.At.tair.db,  # Arabidopsis example
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# KEGG pathway enrichment
ekegg <- enrichKEGG(
  gene = promoter_genes,
  organism = "ath",  # Arabidopsis KEGG code
  pvalueCutoff = 0.05
)
```

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| "TxDb not found" | No annotation database | Create from GFF3 with makeTxDbFromGFF() |
| Many "Intergenic" peaks | tssRegion too narrow | Expand to c(-3000, 1000) |
| Duplicate gene IDs | Multiple transcripts | Use `level="gene"` |
| Slow annotation | Large peak set | Split peaks by chromosome |
| Missing annotations | GFF3 format issues | Validate GFF3 with AGAT |

## Related Tools

| Tool | Purpose |
|------|---------|
| **HOMER annotatePeaks** | Alternative annotation tool |
| **BEDTools closest** | Find nearest genes |
| **ChIPseeker** | R/Bioconductor annotation |
| **GREAT** | Web-based regulatory annotation |

## References

- Yu et al. (2015) ChIPseeker: an R/Bioconductor package for ChIP peak annotation
- Plant promoter architecture: Yamamoto et al. (2009)
