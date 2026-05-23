# bcftools mpileup — Variant Calling for Small Cohorts

**Goal:** SNP/InDel calling via samtools mpileup + bcftools call
**Best for:** Small sample cohorts (<10), non-model species, simple pipelines

## Prerequisites
- samtools 1.10+, bcftools 1.10+
- Sorted, indexed BAM files
- Reference genome with `.fai` index

## Basic Usage

### Standard Pipeline

```bash
# Single-sample calling
bcftools mpileup -f reference.fa sample1_dedup.bam \
  | bcftools call -mv -Oz -o sample1_raw.vcf.gz

# Multi-sample calling
bcftools mpileup -f reference.fa \
  sample1_dedup.bam sample2_dedup.bam sample3_dedup.bam \
  | bcftools call -mv -Oz -o cohort_raw.vcf.gz

# With filtering
bcftools mpileup -f reference.fa \
  -q 20 -Q 20 \            # min mapping quality, min base quality
  sample*.bam \
  | bcftools call -mv -Oz -o raw.vcf.gz

# Index the VCF
bcftools index raw.vcf.gz
```

### Advanced: mpileup + call options

```bash
# Conservative calling (high quality only)
bcftools mpileup -f ref.fa -q 30 -Q 30 --max-depth 10000 sample.bam \
  | bcftools call -mv -P 1e-6 -Oz -o high_conf.vcf.gz

# Polyploid calling (use freebayes for better polyploid support)
bcftools mpileup -f ref.fa sample.bam \
  | bcftools call -mv --ploidy 4 -Oz -o tetraploid.vcf.gz
```

## Key Parameters

| Parameter | Tool | Description | Recommended |
|-----------|------|-------------|-------------|
| `-q` | mpileup | Min mapping quality | 20 (reduce to 10 for polyploids) |
| `-Q` | mpileup | Min base quality | 20 (30 for high-confidence) |
| `--max-depth` | mpileup | Max per-sample depth | 10000 (increase for organelle genomes) |
| `-m` | call | Multiallelic caller | Use for variant sites |
| `-v` | call | Output variant sites only | Use with `-m` |
| `-P` | call | Prior P(ref is not variant) | 1e-6 for diploid, adjust for ploidy |
| `-C 50` | mpileup | Adjust mapping quality for excessive mismatches | Use for high-diversity species |

## When bcftools Over GATK

- **Sample count < 10**: bcftools mpileup is faster and more straightforward than setting up GenomicsDB
- **Non-model species**: No need for known-sites or BQSR training data
- **Limited computational resources**: bcftools uses much less RAM than GATK
- **Quick exploratory analysis**: When you need a first-pass VCF for QC decisions
- **Pooled sequencing**: bcftools supports pooled samples via `-P` settings

## When NOT to Use bcftools

- **Large cohorts (≥30 samples)**: GATK GVCF joint-calling has better sensitivity
- **Low-frequency variant discovery**: GATK's joint calling is more powerful for rare variants
- **When BQSR is needed**: If base quality scores are systematically biased, GATK's BQSR significantly improves results

## Plant-Specific Notes

- **High diversity species**: Plants like maize or sunflower have high SNP density. Use `-C 50` (coefficient for downgrading mapping quality) to avoid overcalling at highly variable regions.
- **Polyploid genotyping**: bcftools mpileup `--ploidy` is functional but less accurate than freebayes or GATK for polyploid genotyping beyond tetraploid level. For autopolyploids, consider freebayes.
- **Chloroplast/mitochondria**: Organelle genomes can have very high coverage (>1000x). Set `--max-depth` high enough to include these but low enough to filter out collapsed repeats. Consider calling organelle variants separately from nuclear variants.
- **Pooled sequencing**: For bulked segregant analysis (BSA) or pool-seq, bcftools mpileup with `-P` adjusted for pool size is the standard approach in plant genomics.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `[E::fai_build] Failed to open` | Missing index | Run `samtools faidx reference.fa` |
| `mpileup` returns empty | No reads overlap | Check BAM coordinate sort and contig names |
| `bcftools call` hangs | High depth at collapsed repeats | Set `--max-depth` |
| VCF has no variants | Threshold too strict | Lower `-P` value or use `-v` flag |
| `[W::vcf_parse] contig not found` | Contig name mismatch | Verify BAM header matches reference |
