# Metagenomics Preprocessing

**Goal:** Quality control, adapter trimming, host read removal for plant metagenomic samples
**Best for:** Preparing raw metagenomic reads for assembly and downstream analysis

## Prerequisites

- fastp (https://github.com/OpenGene/fastp)
- FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Bowtie2 or BWA for host alignment
- Plant reference genome FASTA files

## FastQC — Initial Quality Assessment

```bash
# Run FastQC on all samples
fastqc -t 16 -o fastqc_reports/ raw_reads/*.fastq.gz

# Aggregate with MultiQC
multiqc -o multiqc_report/ fastqc_reports/
```

## fastp — Trimming and Filtering

```bash
# Standard paired-end trimming
fastp -i sample_R1.fq.gz -I sample_R2.fq.gz \
  -o clean_R1.fq.gz -O clean_R2.fq.gz \
  --html fastp_report.html \
  --json fastp_report.json \
  --thread 8 \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --cut_front --cut_tail \
  --cut_window_size 4 --cut_mean_quality 20 \
  --detect_adapter_for_pe

# For low-biomass samples (endosphere/phyllosphere)
fastp -i sample_R1.fq.gz -I sample_R2.fq.gz \
  -o clean_R1.fq.gz -O clean_R2.fq.gz \
  --thread 8 \
  --qualified_quality_phred 30 \
  --length_required 75 \
  --low_complexity_filter \
  --detect_adapter_for_pe
```

## Host Removal with Bowtie2

```bash
# Step 1: Build plant genome index
bowtie2-build plant_genome.fa plant_index

# Step 2: Align reads to plant genome
bowtie2 -x plant_index \
  -1 clean_R1.fq.gz -2 clean_R2.fq.gz \
  -S mapped_to_host.sam \
  --very-sensitive \
  -p 16

# Step 3: Extract unmapped (microbial) reads
samtools view -b -f 12 -F 256 mapped_to_host.sam > microbial.bam

# Step 4: Convert back to FASTQ
samtools fastq -1 microbial_R1.fq -2 microbial_R2.fq microbial.bam
gzip microbial_R1.fq microbial_R2.fq
```

## KneadData (Integrated QC + Host Removal)

```bash
# Create host database
kneaddata_database --download human_genome bowtie2 db/

# For plant host: build custom database
bowtie2-build plant_genome.fa plant_db/plant_genome

# Run KneadData
kneaddata -i sample_R1.fq.gz -i sample_R2.fq.gz \
  -o kneaddata_output/ \
  -db plant_db/plant_genome \
  --trimmomatic /path/to/Trimmomatic \
  --trimmomatic-options \
    "SLIDINGWINDOW:4:20 MINLEN:50" \
  -t 16 \
  --run-fastqc-start --run-fastqc-end
```

## Quality Report

```bash
# Generate comprehensive QC report
python3 << 'EOF'
import json, os

for fq in ["clean_R1.fq.gz"]:  # per sample
    with open("fastp_report.json") as f:
        stats = json.load(f)
    print(f"Total reads: {stats['summary']['before_filtering']['total_reads']}")
    print(f"After QC: {stats['summary']['after_filtering']['total_reads']}")
    print(f"Retained: {stats['summary']['after_filtering']['total_reads'] / stats['summary']['before_filtering']['total_reads'] * 100:.1f}%")

# Host removal check
total_microbial = sum(1 for _ in open("microbial_R1.fq")) / 4
print(f"Microbial reads: {total_microbial}")
EOF
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `--qualified_quality_phred` | 20 (soil) / 30 (endosphere) | Stricter for low-biomass |
| `--length_required` | 50 (soil) / 75 (endosphere) | Longer for low-biomass |
| Bowtie2 `--very-sensitive` | Enabled | Maximize host removal |
| `--low_complexity_filter` | Endosphere only | Reduce sequencing artifacts |

## Plant-Specific Notes

- Use combined reference: nuclear + chloroplast + mitochondrial genome for host removal
- If plant reference is incomplete: supplement with Kraken2 classification-based filtering
- Check host removal efficiency: `grep -c "^+" microbial_R1.fq` vs original
- Rhizosphere: expect 10-30% plant reads; Endosphere: 80-95% (before removal)

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Zero microbial reads after host removal | All reads mapped to host | Check if sample is pure plant tissue |
| Adapter not detected | Custom library prep | Manually specify adapter sequences |
| Bowtie2 index too large | Plant genome >3 Gb | Use BWA-MEM instead, or split chromosomes |
| fastp slow | Too many threads for I/O | Reduce `-t` to match disk I/O |
