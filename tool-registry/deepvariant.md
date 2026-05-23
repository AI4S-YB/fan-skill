# DeepVariant — Deep Learning Variant Caller

**Goal:** SNP/InDel calling using a convolutional neural network
**Best for:** Long-read data (PacBio HiFi, ONT), when maximum accuracy is needed

## Prerequisites
- DeepVariant 1.4+ (Docker or pre-built binary)
- Sorted, indexed BAM file
- Reference genome with `.fai`
- GPU recommended but not required (CPU mode works, 10-20x slower)

## Basic Usage

```bash
# Using Docker (recommended for easy setup)
BIN_VERSION="1.5.0"

# PacBio HiFi
docker run -v /path/to/data:/data \
  google/deepvariant:${BIN_VERSION} \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref /data/reference.fa \
  --reads /data/sample1_sorted.bam \
  --output_vcf /data/sample1_deepvariant.vcf.gz \
  --output_gvcf /data/sample1_deepvariant.g.vcf.gz \
  --num_shards 8

# ONT
docker run -v /path/to/data:/data \
  google/deepvariant:${BIN_VERSION} \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type ONT_R104 \
  --ref /data/reference.fa \
  --reads /data/sample1_sorted.bam \
  --output_vcf /data/sample1_deepvariant.vcf.gz \
  --num_shards 8

# Illumina (also supported)
docker run -v /path/to/data:/data \
  google/deepvariant:${BIN_VERSION} \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WGS \
  --ref /data/reference.fa \
  --reads /data/sample1_sorted.bam \
  --output_vcf /data/sample1_deepvariant.vcf.gz \
  --num_shards 8
```

## Key Parameters

| Parameter | Description | Recommendation |
|-----------|-------------|----------------|
| `--model_type` | Sequencing technology | PACBIO / ONT_R104 / WGS / WES |
| `--num_shards` | Parallel shards for speed | Match CPU core count |
| `--output_gvcf` | gVCF output | Use for multi-sample joint calling |
| `--regions` | Limit to specific contigs | Use for large genomes to split work |
| `--vsc_min_fraction` | Variant support confidence | 0.12 (default); increase for stringency |

## When DeepVariant

- **PacBio HiFi data**: DeepVariant achieves consistently higher accuracy than GATK for PacBio HiFi reads, particularly for InDels
- **ONT data**: The ONT_R104 model is specifically trained on nanopore error profiles
- **Single-sample analysis**: DeepVariant works well per-sample; multi-sample can be done with GLnexus
- **When you need maximum accuracy**: DeepVariant's neural network model picks up subtle signal patterns missed by statistical callers

## When NOT DeepVariant

- **Large sample cohorts without GPU**: DeepVariant is 10-20x slower on CPU than GATK HaplotypeCaller
- **Polyploid genotyping**: DeepVariant currently only supports diploid calling; cannot handle polyploid dosage
- **Non-model species with no training data**: DeepVariant training data is human/mammalian-heavy; accuracy on highly divergent species may be slightly reduced (though still good for PacBio HiFi)
- **Limited computational resources**: GATK or bcftools are much lighter

## Plant-Specific Notes

- **Applicability to plant genomes**: DeepVariant works well on plant PacBio HiFi data. The model generalizes reasonably well across species because the image-based pileup representation is species-agnostic. However, for species with extreme GC content or very large repetitive genomes, benchmark against GATK before committing.
- **Polyploid limitation**: DeepVariant does not support polyploid calling (>2 copies). For wheat, potato, or other polyploid crops, use GATK with `-ploidy` or freebayes.
- **Organelle variants**: DeepVariant may overcall variants in high-coverage organelle regions. Filter organelle contigs and call separately with appropriate depth expectations.
- **GLnexus for multi-sample**: Use GLnexus (https://github.com/dnanexus-rnd/GLnexus) to merge per-sample DeepVariant gVCFs into a joint-called VCF, analogous to GATK GenotypeGVCFs.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `docker: command not found` | Docker not installed | Use pre-built binary or Singularity |
| `CUDA not available` | No GPU in container | Uses CPU automatically (slower but works) |
| `Could not find reference index` | Missing .fai | Run `samtools faidx reference.fa` |
| `Out of memory` | Genome too large for shard size | Reduce `--num_shards` or use `--regions` |
| `make_examples` phase slow | Large repetitive regions | Limit to non-repetitive regions for first pass |
