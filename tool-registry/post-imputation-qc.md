# Post-Imputation Quality Control

**Goal:** Filter imputed variants by quality metrics, retaining high-confidence calls for downstream analysis
**Best for:** Output from Beagle, IMPUTE2, or Minimac4 imputation

## Prerequisites

- Imputed VCF with embedded quality metrics (DR2, INFO, R²)
- PLINK 1.9/2.0 or bcftools
- Decision on downstream analysis type (GWAS, GS, QTL mapping, etc.)

## Standard Post-Imputation Filtering

### For Beagle Output

Beagle stores dosage R-squared (DR2) in the VCF INFO field.

```bash
# Step 1: Filter by DR2 >= 0.3
bcftools filter -i 'DR2>=0.3' imputed.vcf.gz -Oz -o step1_dr2.vcf.gz

# Step 2: Filter by MAF >= 0.01
bcftools filter -i 'AF>=0.01 & AF<=0.99' step1_dr2.vcf.gz -Oz -o step2_maf.vcf.gz

# Step 3: Convert to PLINK for further QC
plink2 --vcf step2_maf.vcf.gz \
  --make-bed --out imputed_filtered

# Step 4: Check missingness in imputed data (should be very low)
plink2 --bfile imputed_filtered \
  --geno 0.05 \
  --maf 0.01 \
  --make-bed --out imputed_final
```

### For IMPUTE2 Output

IMPUTE2 produces a separate info file. Filter variants with info < 0.3.

```bash
# Step 1: Extract variants with info >= 0.3
awk '$5 >= 0.3' output_info > good_variants.txt

# Step 2: Keep only good variants
plink2 --gen output.gen \
  --sample samples.txt \
  --extract good_variants.txt \
  --maf 0.01 \
  --make-bed --out imputed_filtered
```

### For Minimac4 Output

Minimac4 includes R² in the info file or VCF INFO field.

```bash
# Step 1: Filter by R² >= 0.3
bcftools filter -i 'R2>=0.3' imputed.dose.vcf.gz -Oz -o step1_r2.vcf.gz

# Step 2: Filter by MAF
bcftools filter -i 'MAF>=0.01' step1_r2.vcf.gz -Oz -o imputed_filtered.vcf.gz

# Step 3: Convert to PLINK
plink2 --vcf imputed_filtered.vcf.gz \
  --make-bed --out imputed_final
```

## Strict Filtering (for GWAS downstream)

When downstream analysis is GWAS, apply stricter thresholds:

```bash
# GWAS-strict filtering
bcftools filter -i 'DR2>=0.5 && AF>=0.05 && AF<=0.95' imputed.vcf.gz -Oz -o imputed_strict.vcf.gz

# Alternative: using PLINK2 with hard-call threshold
plink2 --vcf imputed.vcf.gz \
  --hard-call-threshold 0.9 \  # Convert dosage to hard calls with high confidence
  --maf 0.05 \
  --geno 0.02 \
  --hwe 1e-6 \
  --make-bed --out imputed_gwas_ready
```

## Genomic Selection (GS) Post-Imputation

For GS downstream, retain more markers even at lower quality:

```bash
# GS-relaxed filtering
bcftools filter -i 'DR2>=0.3 && AF>=0.01' imputed.vcf.gz -Oz -o imputed_gs.vcf.gz

# Missing genotypes can be mean-imputed in GS models
# but prefer to have DR2 >= 0.3 as a minimum
```

## Quality Metrics to Report

After post-imputation QC, generate these metrics:

| Metric | Command | Good | Acceptable | Poor |
|--------|---------|------|------------|------|
| Mean R²/DR2 | bcftools query | > 0.8 | 0.6-0.8 | < 0.6 |
| % variants R² > 0.3 | Custom script | > 90% | 80-90% | < 80% |
| % variants R² > 0.8 | Custom script | > 60% | 40-60% | < 40% |
| Post-QC variant count | bcftools index --nrecords | > input SN | 50-100% input | < 50% input |
| Post-QC sample count | bcftools query -l | = input N | = input N | < input N |

## R² by MAF Bin

Imputation accuracy varies by allele frequency. Check per-bin:

```bash
# Extract DR2 and AF for all variants
bcftools query -f '%CHROM\t%POS\t%AF\t%DR2\n' imputed.vcf.gz > variant_metrics.txt

# Use R to summarize by MAF bin
Rscript -e '
dat <- read.table("variant_metrics.txt", header=FALSE,
                   col.names=c("CHR","POS","AF","DR2"))
dat$MAF <- pmin(dat$AF, 1-dat$AF)
dat$MAF_bin <- cut(dat$MAF, breaks=c(0, 0.01, 0.05, 0.1, 0.2, 0.5))
aggregate(DR2 ~ MAF_bin, dat, function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
'
```

**Interpretation**: If DR2 is very low (<0.3) for MAF < 0.01, this is expected and normal. Flag only if DR2 is consistently low across all MAF bins.

## Per-Chromosome QC

Check for chromosome-level issues:

```bash
# Mean DR2 per chromosome
for chr in $(seq 1 12); do
  bcftools view -r "Chr${chr}" imputed.vcf.gz | \
    bcftools query -f '%DR2\n' | \
    awk -v chr="$chr" '{sum+=$1; n++} END{print chr, sum/n}'
done
```

**Red flag**: A single chromosome with substantially lower mean R² than others. Possible cause: poor reference panel coverage for that chromosome, or a misassembled region.

## Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| Most variants DR2 < 0.3 | Failed imputation | Check input quality; re-run with different tool |
| Certain MAF bin has low DR2 | Reference panel coverage gap | Expected for very rare / very common variants |
| One chromosome outlier | Bad reference or assembly | Investigate that chromosome; may need to exclude |
| Too few variants retained under strict filter | Standard thresholds too tight | Use moderate thresholds if GWAS is exploratory |
| bcftools filter fails | DR2 field name varies by tool | Check VCF header: DR2 (Beagle), INFO (IMPUTE2), R2 (Minimac4) |
