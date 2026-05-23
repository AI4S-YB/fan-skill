# BWA-MEM & minimap2 — Sequence Alignment

**Goal:** Map resequencing reads to a reference genome
**Best for:** BWA-MEM for Illumina short reads; minimap2 for PacBio/ONT long reads

## Prerequisites
- Reference genome with index (`.fa` + `.fai` + `.dict`)
- For BWA-MEM: `bwa index reference.fa`
- For minimap2: no separate index step needed
- Quality-trimmed FASTQ files

## Basic Usage

### BWA-MEM (Illumina short reads)

```bash
# Index reference (one-time)
bwa index reference.fa
samtools faidx reference.fa
gatk CreateSequenceDictionary -R reference.fa

# Align paired-end reads
bwa mem -t 8 -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1" \
  reference.fa \
  sample1_R1.fastq.gz sample1_R2.fastq.gz \
  | samtools sort -@ 4 -o sample1_sorted.bam -

samtools index sample1_sorted.bam
```

### minimap2 (PacBio/ONT long reads)

```bash
# Align PacBio HiFi reads
minimap2 -t 8 -ax map-hifi reference.fa sample1.fastq.gz \
  | samtools sort -@ 4 -o sample1_sorted.bam -

# Align ONT reads
minimap2 -t 8 -ax map-ont reference.fa sample1.fastq.gz \
  | samtools sort -@ 4 -o sample1_sorted.bam -

samtools index sample1_sorted.bam
```

## Key Parameters

| Parameter | Tool | Description | Recommended |
|-----------|------|-------------|-------------|
| `-t` | bwa/minimap2 | Threads | Match available cores |
| `-R` | bwa mem | Read group | Always include; required by GATK |
| `-M` | bwa mem | Mark shorter split hits as secondary | Use for compatibility |
| `-x map-hifi` | minimap2 | PacBio HiFi preset | Use for CCS/HiFi reads |
| `-x map-ont` | minimap2 | ONT preset | Use for nanopore reads |
| `-ax sr` | minimap2 | Short-read preset | Alternative for Illumina |

## QC After Alignment

```bash
# Mapping statistics
samtools flagstat sample1_sorted.bam
samtools stats sample1_sorted.bam | grep ^SN

# Coverage
samtools depth -a sample1_sorted.bam | awk '{sum+=$3; n++} END {print sum/n}'
```

## Mark Duplicates

```bash
# Using GATK (required for GATK pipeline)
gatk MarkDuplicates \
  -I sample1_sorted.bam \
  -O sample1_dedup.bam \
  -M sample1_metrics.txt

# Using samtools ( lighter, adequate for bcftools pipeline)
samtools markdup -r sample1_sorted.bam sample1_dedup.bam
```

## Plant-Specific Notes

- **Polyploid mapping rates**: Expect lower mapping rates in polyploids (70-85%) compared to diploids (>90%). Reads mapping to homeologous regions may have low MAPQ — do not be alarmed by a 75% mapping rate in hexaploid wheat.
- **Organelle genomes**: Chloroplast and mitochondrial reads can constitute 5-20% of total reads. If interested in organelle variants, align separately to organelle references; otherwise filter out organelle-mapped reads to avoid inflated coverage in those regions.
- **Large genomes**: Some plant genomes (wheat ~17Gb, pine ~22Gb) require substantial RAM. Consider `bwa mem -c` to cap the number of hits, or split chromosomes for parallel alignment.
- **Read groups**: Always tag with proper `@RG`; GATK HaplotypeCaller requires it. The `SM` (sample name) tag is essential for multi-sample calling.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `[E::bwa_idx_load] fail to locate the index files` | Missing BWA index | Run `bwa index reference.fa` |
| `Read group information is missing` | No @RG tag | Add `-R` to bwa mem command |
| `MATE_NOT_FOUND` warnings | Discordant read pairs in polyploid | Normal for polyploids; filter with `samtools view -f 2` if needed |
| Low mapping rate (<60%) | Wrong reference or contamination | Check species identity; run FastQC |
| `[mem_sam_pe] paired reads have different names` | FASTQ order mismatch | Check that R1/R2 files are properly paired |
