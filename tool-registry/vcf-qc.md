# VCF QC & Filtering — Hard Filters and VQSR

**Goal:** Filter raw variant calls to a high-confidence set
**Best for:** All VCF pipelines; choose hard filtering or VQSR based on sample count

## Prerequisites
- Raw VCF from GATK/bcftools/DeepVariant
- bcftools or GATK
- Known sites VCF (for VQSR only)

## Hard Filtering

### GATK Best Practices Hard Filters (for SNP)

```
QD < 2.0          → low quality by depth
FS > 60.0         → strand bias (Fisher strand)
MQ < 40.0         → low mapping quality
MQRankSum < -12.5 → mapping quality rank sum
ReadPosRankSum < -8.0 → read position rank sum
SOR > 3.0         → strand odds ratio
```

### GATK Best Practices Hard Filters (for INDEL)

```
QD < 2.0
FS > 200.0
ReadPosRankSum < -20.0
SOR > 10.0
```

### Applying Hard Filters with GATK

```bash
# Extract SNPs
gatk SelectVariants \
  -R reference.fa \
  -V raw.vcf.gz \
  --select-type-to-include SNP \
  -O raw_snps.vcf.gz

# Apply SNP filters
gatk VariantFiltration \
  -R reference.fa \
  -V raw_snps.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  -O filtered_snps.vcf.gz

# Extract INDELs
gatk SelectVariants \
  -R reference.fa \
  -V raw.vcf.gz \
  --select-type-to-include INDEL \
  -O raw_indels.vcf.gz

# Apply INDEL filters
gatk VariantFiltration \
  -R reference.fa \
  -V raw_indels.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  --filter-expression "SOR > 10.0" --filter-name "SOR10" \
  -O filtered_indels.vcf.gz

# Merge and keep only PASS
gatk MergeVcfs \
  -I filtered_snps.vcf.gz \
  -I filtered_indels.vcf.gz \
  -O combined_filtered.vcf.gz

bcftools view -f PASS combined_filtered.vcf.gz -Oz -o final_pass.vcf.gz
```

### Hard Filtering with bcftools

```bash
bcftools filter -i 'QD > 2.0 && FS < 60 && MQ > 40 && SOR < 3.0' \
  -Oz -o filtered.vcf.gz raw.vcf.gz
```

## VQSR (Variant Quality Score Recalibration)

**Requirement**: ≥30 samples + known sites VCF

```bash
# Build SNP recalibration model
gatk VariantRecalibrator \
  -R reference.fa \
  -V raw.vcf.gz \
  --resource:known,known=false,training=true,truth=true,prior=10 known_sites.vcf.gz \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
  -mode SNP \
  -O snp_recal.table \
  --tranches-file snp_tranches

# Apply SNP VQSR
gatk ApplyVQSR \
  -R reference.fa \
  -V raw.vcf.gz \
  -O recalibrated_snps.vcf.gz \
  --truth-sensitivity-filter-level 99.7 \
  --tranches-file snp_tranches \
  --recal-file snp_recal.table \
  -mode SNP

# Repeat for INDEL with -mode INDEL
```

## When VQSR vs Hard Filtering

| Criterion | Hard Filtering | VQSR |
|-----------|---------------|------|
| Sample count | Any | ≥ 30 |
| Known sites required | No | Yes |
| Accuracy | Good | Better (at scale) |
| Simplicity | Simple, transparent | Complex, black-box |
| False positive control | OK | Excellent |
| False negative risk | Higher | Lower |

## QC Metrics After Filtering

```bash
# Ti/Tv ratio (expect 2.0-2.5 for plants)
bcftools stats final_pass.vcf.gz | grep "TSTV"

# Missing rate per sample
bcftools stats -s - final_pass.vcf.gz | grep "PSC"

# Heterozygosity per sample
bcftools stats -s - final_pass.vcf.gz | grep "PSC"

# Per-sample depth
bcftools stats -s - final_pass.vcf.gz | grep "DP"
```

## Plant-Specific Notes

- **Ti/Tv ratio expectations**: Plants have lower Ti/Tv than humans (expected ~2.0-2.5 in coding regions, 1.5-2.0 genome-wide). Values below 1.5 suggest excessive false positive calls. Values above 3.0 suggest missed variants.
- **Polyploid filtering**: Standard hard filter thresholds assume diploid data. For polyploids, expect higher FS values (reads from homeologous regions cause apparent strand bias). Relax FS threshold by ~20-30% for polyploids.
- **Inbred species**: Very low heterozygosity is expected. Do not filter on heterozygosity for inbred crops like rice, soybean, or wheat. A near-zero heterozygosity rate is a sign of data quality, not a problem.
- **Missing data in polyploids**: Expect higher missing rates in polyploids because reads that map equally well to multiple subgenomes result in low-confidence calls. Do not over-filter on missingness for polyploid datasets.
- **Relaxed thresholds for low-coverage data**: For 5-10x coverage, relax QD threshold to 1.5 and MQ to 30. Low-coverage plant resequencing is common for population surveys.
- **Organelle variant filtering**: Chloroplast and mitochondrial variants can confound Ti/Tv and heterozygosity metrics. Filter organelle contigs out before computing genome-wide QC stats, or compute stats separately per compartment.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| VQSR fails to build model | Too few variants or samples | Fall back to hard filtering |
| `IllegalArgumentException: Badly formed -an` | Incorrect annotation name | Check available annotations: `gatk VariantRecalibrator --list` |
| Filtered out 99% of variants | Threshold too strict | Check raw data quality; relax filter thresholds |
| Ti/Tv < 1.0 | Massive false positives or sequencing error | Re-examine raw data quality and alignment; check for contamination |
| bcftools filter: no variants pass | Expression syntax error | Use `-i` for include, `-e` for exclude; test with `--dry-run` |
