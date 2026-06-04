# Quality Control: ChIP-seq/ATAC-seq Quality Assessment

## Tool Overview

**Tool ID**: `quality-control`
**Category**: Epigenomics Quality Control
**Purpose**: Assess and validate ChIP-seq/ATAC-seq data quality before downstream analysis
**Primary Tools**: FastQC, Picard, deepTools, ChIPQC
**Alternative**: NF-Core QC modules, MultiQC

## When to Use

- Pre-alignment quality assessment (FastQC)
- Post-alignment QC (duplication, enrichment)
- Peak-level quality metrics (FRiP, NSC, RSC)
- Sample comparison and batch effect detection

## Key Parameters

### ChIP-seq Quality Metrics

| Metric | Threshold (Good) | Threshold (Acceptable) | Description |
|--------|------------------|------------------------|-------------|
| **FRiP** (Fraction of Reads in Peaks) | > 0.10 (TF) | > 0.05 (TF) | Fraction of reads overlapping peaks |
| **FRiP** (Histone) | > 0.30 | > 0.15 | Higher for broad marks |
| **NSC** (Normalized Strand Cross-correlation) | > 1.10 | > 1.05 | Peak enrichment signal |
| **RSC** (Relative Strand Cross-correlation) | > 0.80 | > 0.50 | Library complexity indicator |
| **NRF** (Non-Redundant Fraction) | > 0.80 | > 0.50 | Unique reads / total reads |
| **PBC** (PCR Bottleneck Coefficient) | > 0.80 | > 0.50 | Library complexity measure |

### ATAC-seq Specific Metrics

| Metric | Threshold (Good) | Threshold (Acceptable) | Description |
|--------|------------------|------------------------|-------------|
| **TSS Enrichment** | > 10 | > 5 | Enrichment at transcription start sites |
| **Nucleosome-free reads** | > 50% | > 30% | Fraction of <100bp fragments |
| **Organelle contamination** | < 5% | < 15% | Chloroplast/mitochondrial reads |
| **Fragment size periodicity** | Clear 200bp | Visible 200bp | Nucleosome pattern |

## Recommended Parameter Sets

### ChIP-seq QC with ChIPQC

```r
library(ChIPQC)

# Create sample sheet
samples <- data.frame(
  sampleID = c("WT_1", "WT_2", "KO_1", "KO_2"),
  tissue = "leaf",
  factor = "TF",
  condition = c("WT", "WT", "KO", "KO"),
  replicate = c(1, 2, 1, 2),
  bamReads = c("WT_1.bam", "WT_2.bam", "KO_1.bam", "KO_2.bam"),
  Peaks = c("WT_1_peaks.narrowPeak", "WT_2_peaks.narrowPeak", 
            "KO_1_peaks.narrowPeak", "KO_2_peaks.narrowPeak")
)

# Run QC
chipqc <- ChIPQC(samples, consensus = TRUE, bCount = TRUE)

# Generate report
ChIPQCreport(chipqc, reportName = "ChIPQC_report")
```

### ATAC-seq QC with deepTools

```bash
# 1. TSS enrichment calculation
computeMatrix reference-point \
  -S sample.bigWig \
  -R genes.bed \
  --referencePoint TSS \
  -b 2000 -a 2000 \
  -out matrix.gz

plotProfile -m matrix.gz -out TSS_profile.png

# 2. Fragment size distribution
plotFingerprint \
  -b sample.bam \
  -plot fingerprint.png \
  --outRawCounts counts.txt
```

## Plant-Specific Quality Considerations

### 1. Chloroplast/Mitochondrial Contamination

Plant samples often have organelle contamination:

```bash
# Check organelle contamination
samtools idxstats sample.bam | grep -E "chloroplast|mitochondria" | \
  awk '{sum+=$3} END {print "Organelle reads: " sum}'

# Filter out organelle reads (recommended for ATAC-seq)
samtools view -b -L nuclear_chromosomes.bed sample.bam > sample_filtered.bam
```

**Typical plant organelle contamination levels:**
- Good: < 5% organelle reads
- Acceptable: 5-15% organelle reads
- Problematic: > 15% organelle reads

### 2. FRiP Expectations for Plants

Plant ChIP-seq typically has lower FRiP than animal studies:

| Sample Type | Expected FRiP | Reason |
|-------------|---------------|--------|
| **Plant TF** | 0.05 - 0.15 | Cell wall limits crosslinking efficiency |
| **Plant Histone** | 0.15 - 0.35 | Broad marks easier to capture |
| **Animal TF** | 0.10 - 0.25 | Higher crosslinking efficiency |
| **ATAC-seq** | 0.20 - 0.50 | Depends on tissue and nuclei isolation |

### 3. TSS Enrichment for Plant ATAC-seq

```r
# Calculate TSS enrichment (plant-specific thresholds)
calculate_tss_enrichment <- function(bam_file, tss_bed) {
  # Get reads at TSS ± 2000bp
  tss_reads <- count_reads_in_region(bam_file, tss_bed, extend = 2000)
  
  # Get reads in flanking region (± 2000-3000bp)
  flank_reads <- count_reads_in_region(bam_file, tss_bed, 
                                         inner = 2000, outer = 3000)
  
  # Calculate enrichment
  enrichment <- mean(tss_reads) / mean(flank_reads)
  
  # Plant thresholds
  if (enrichment > 10) return("Excellent")
  if (enrichment > 5) return("Good")
  if (enrichment > 2) return("Acceptable")
  return("Poor")
}
```

## QC Pipeline Integration

### Complete Pre-alignment QC

```bash
# FastQC on raw reads
fastqc -t 8 sample_R1.fastq.gz sample_R2.fastq.gz

# MultiQC summary
multiqc . -o multiqc_report

# Check for adapters
# Look for "Adapter Content" module in FastQC report
```

### Post-alignment QC Pipeline

```bash
# 1. Alignment statistics
samtools flagstat sample.bam > flagstat.txt

# 2. Duplication rate (Picard)
picard MarkDuplicates \
  I=sample.bam \
  O=dedup.bam \
  M=dup_metrics.txt \
  REMOVE_DUPLICATES=true

# 3. Insert size distribution
picard CollectInsertSizeMetrics \
  I=sample.bam \
  O=insert_metrics.txt \
  H=insert_histogram.pdf

# 4. Coverage uniformity
samtools depth sample.bam | \
  awk '{sum+=$3; sumsq+=$3*$3} END {print "Mean depth: " sum/NR}'

# 5. Peak calling for FRiP calculation
macs2 callpeak -t sample.bam -f BAM -g 1.2e8 -n sample --outdir peaks/

# 6. Calculate FRiP
total_reads=$(samtools view -c sample.bam)
peak_reads=$(bedtools intersect -u -a sample.bam -b peaks/sample_peaks.narrowPeak | samtools view -c -)
frip=$(echo "scale=4; $peak_reads / $total_reads" | bc)
echo "FRiP: $frip"
```

## Quality Control Report Template

```r
# Generate comprehensive QC report
generate_qc_report <- function(sample_info) {
  
  report <- list(
    # Alignment metrics
    alignment_rate = get_alignment_rate(sample_info$bam),
    duplication_rate = get_duplication_rate(sample_info$bam),
    
    # Enrichment metrics
    frip = calculate_frip(sample_info$bam, sample_info$peaks),
    nsc_rsc = calculate_cross_correlation(sample_info$bam),
    
    # Library complexity
    nrf = calculate_nrf(sample_info$bam),
    pbc = calculate_pbc(sample_info$bam),
    
    # Sample-specific
    tss_enrichment = if(sample_info$type == "ATAC") 
      calculate_tss_enrichment(sample_info$bam),
    organelle_contamination = if(sample_info$species == "plant")
      calculate_organelle_contamination(sample_info$bam)
  )
  
  # Quality flags
  report$quality_flag <- assess_quality(report)
  
  return(report)
}
```

## Troubleshooting Common Issues

### Low FRiP

| Problem | Possible Causes | Solutions |
|---------|-----------------|-----------|
| FRiP < 0.01 | Poor IP efficiency | Check antibody quality; optimize IP conditions |
| FRiP = 0.01-0.05 | Weak enrichment | Increase sequencing depth; check antibody |
| FRiP varies between replicates | Inconsistent protocol | Standardize IP conditions |

### Low Library Complexity

| NRF Value | Interpretation | Action |
|-----------|----------------|--------|
| NRF < 0.3 | Severe PCR bottleneck | Reduce PCR cycles; re-sequence with more input |
| NRF 0.3-0.5 | Moderate bottleneck | Consider deduplication impact on analysis |
| NRF > 0.5 | Good complexity | Proceed with analysis |

### High Organelle Contamination (ATAC-seq)

```bash
# Solutions for high organelle contamination:
# 1. Improve nuclei isolation
# 2. Add chloroplast lysis step
# 3. Use chloroplast-depleted tissue (e.g., root vs leaf)
# 4. Filter bioinformatically (last resort)

# Filter organelle reads
samtools view -h sample.bam | \
  grep -v "chloroplast\|mitochondria" | \
  samtools view -b - > sample_nuclear_only.bam
```

## MultiQC Integration

```bash
# Collect all QC metrics in one report
multiqc . \
  --module fastqc \
  --module picard \
  --module macs2 \
  --module samtools \
  -o multiqc_report
```

## Quality Thresholds Summary

### For ChIP-seq (TF)

| Metric | Pass | Warn | Fail |
|--------|------|------|------|
| FRiP | > 0.05 | 0.01-0.05 | < 0.01 |
| NSC | > 1.10 | 1.05-1.10 | < 1.05 |
| RSC | > 0.80 | 0.50-0.80 | < 0.50 |
| NRF | > 0.50 | 0.30-0.50 | < 0.30 |
| Duplicates | < 0.50 | 0.50-0.70 | > 0.70 |

### For ATAC-seq

| Metric | Pass | Warn | Fail |
|--------|------|------|------|
| TSS Enrichment | > 8 | 5-8 | < 5 |
| Nucleosome-free | > 40% | 25-40% | < 25% |
| Organelle % | < 10% | 10-20% | > 20% |
| NRF | > 0.50 | 0.30-0.50 | < 0.30 |

## References

- Landt et al. (2012) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia
- ENCODE Quality Metrics: https://www.encodeproject.org/chip-seq/transcription_factor/
- Plant ATAC-seq: Mahé et al. (2021)
