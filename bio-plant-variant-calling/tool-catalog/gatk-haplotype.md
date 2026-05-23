# GATK4 HaplotypeCaller — Germline Variant Calling

**Goal:** SNP/InDel discovery via HaplotypeCaller + joint genotyping
**Best for:** Diploid species with ≥10 samples, well-annotated reference genome

## Prerequisites
- GATK 4.2+
- Sorted, duplicate-marked BAM files with read groups
- Reference genome with `.dict` and `.fai`
- Known sites VCF (optional, for BQSR)

## Pipeline Overview

```
BAM → BaseRecalibrator → ApplyBQSR → HaplotypeCaller (per-sample GVCF)
  → GenomicsDBImport → GenotypeGVCFs → VCF
```

## Basic Usage

### Step 1: Base Quality Score Recalibration (BQSR)

```bash
# Generate recalibration table
gatk BaseRecalibrator \
  -R reference.fa \
  -I sample1_dedup.bam \
  --known-sites known_sites.vcf.gz \
  -O sample1_recal.table

# Apply recalibration
gatk ApplyBQSR \
  -R reference.fa \
  -I sample1_dedup.bam \
  --bqsr-recal-file sample1_recal.table \
  -O sample1_recal.bam
```

### Step 2: HaplotypeCaller (per-sample GVCF)

```bash
# Diploid
gatk HaplotypeCaller \
  -R reference.fa \
  -I sample1_recal.bam \
  -ERC GVCF \
  -O sample1.g.vcf.gz

# Polyploid (e.g., tetraploid wheat)
gatk HaplotypeCaller \
  -R reference.fa \
  -I sample1_recal.bam \
  -ERC GVCF \
  -ploidy 4 \
  -O sample1.g.vcf.gz
```

### Step 3: Joint Genotyping

```bash
# Import GVCFs into GenomicsDB
gatk GenomicsDBImport \
  -V sample1.g.vcf.gz \
  -V sample2.g.vcf.gz \
  -V sample3.g.vcf.gz \
  --genomicsdb-workspace-path gendb_workspace \
  --intervals chr1.bed

# Joint genotyping
gatk GenotypeGVCFs \
  -R reference.fa \
  -V gendb://gendb_workspace \
  -O cohort_raw.vcf.gz
```

## Key Parameters

| Parameter | Tool | Description | Plant Recommendation |
|-----------|------|-------------|---------------------|
| `-ERC GVCF` | HaplotypeCaller | Emit reference confidence | Always use for multi-sample |
| `-ploidy` | HaplotypeCaller | Ploidy level | Critical for polyploids |
| `--heterozygosity` | HaplotypeCaller | Expected heterozygosity | 0.001 for inbred, 0.01 for outcross |
| `--min-base-quality-score` | HaplotypeCaller | Minimum base quality | 20 (raise to 30 for low-coverage) |
| `--intervals` | GenomicsDBImport | Scaffold-by-scaffold import | Essential for large plant genomes |
| `--max-alternate-alleles` | GenotypeGVCFs | Max alt alleles per site | 6 (increase for polyploid) |

## Plant-Specific Notes

- **Ploidy**: Set `-ploidy` correctly. Default is 2 (diploid). For hexaploid wheat use 6, for tetraploid potato use 4. Incorrect ploidy will produce severely biased genotype calls.
- **Large genomes**: Import GVCFs scaffold-by-scaffold using `--intervals`. Attempting to import all scaffolds into a single GenomicsDB workspace for a large plant genome will exhaust memory.
- **BQSR without known sites**: For non-model plants without known variant databases, run an initial round of variant calling with lenient filters to create a bootstrap set of known sites, then re-run BQSR.
- **Heterozygosity expectations**: Self-pollinated crops (rice, soybean) expect very low heterozygosity (~0.001). Cross-pollinated crops (maize, sunflower) expect higher rates (~0.01). Setting `--heterozygosity` correctly improves indel calling accuracy.
- **Joint calling vs per-sample**: For small cohorts (<10 samples), per-sample calling with `bcftools merge` is simpler and gives comparable results. For large cohorts, GVCF joint-calling gives better sensitivity for rare variants.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `A USER ERROR has occurred: Bad input: Read group is missing` | BAM missing @RG | Re-align with `-R` read group flag |
| `java.lang.OutOfMemoryError` | Java heap too small | Increase with `--java-options "-Xmx32G"` |
| `The interval ... does not exist in the sequence dictionary` | Contig name mismatch | Check reference vs BAM contig names |
| `IllegalArgumentException: ploidy must be > 0` | Wrong syntax | Use `-ploidy 4` not `--ploidy 4` |
| GVCF import very slow | Too many intervals at once | Split by chromosome/scaffold |
