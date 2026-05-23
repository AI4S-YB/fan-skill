# Pre-Imputation Quality Control

**Goal:** Prepare genotype data for imputation by removing low-quality variants and samples
**Best for:** All genotype data before any imputation tool

## Prerequisites

- PLINK 1.9 or PLINK2
- Input VCF or PLINK binary (bed/bim/fam)
- Reference genome info (for strand/allele checking)

## Standard QC Pipeline (snp_count >= 1000)

```bash
# Step 1: Filter by missingness per variant
# Remove variants with >10% missing calls
plink2 --bfile input \
  --geno 0.10 \
  --make-bed --out step1_geno

# Step 2: Filter by missingness per sample
# Remove samples with >10% missing calls
plink2 --bfile step1_geno \
  --mind 0.10 \
  --make-bed --out step2_mind

# Step 3: Filter by MAF
# Remove very rare variants (MAF < 0.01)
plink2 --bfile step2_mind \
  --maf 0.01 \
  --make-bed --out step3_maf

# Step 4: HWE filter (skip for inbred species)
# Remove variants with extreme HWE deviation
plink2 --bfile step3_maf \
  --hwe 1e-6 \
  --make-bed --out step4_hwe

# Step 5: Final conversion to VCF for imputation
plink2 --bfile step4_hwe \
  --recode vcf bgz \
  --out cleaned_for_imputation
```

## Relaxed QC Pipeline (snp_count < 1000)

```bash
# For low-density markers: keep more variants
plink2 --bfile input \
  --geno 0.20 \      # Relaxed: allow 20% missing
  --mind 0.10 \
  --maf 0.01 \
  --make-bed --out cleaned_relaxed

plink2 --bfile cleaned_relaxed \
  --recode vcf bgz \
  --out cleaned_for_imputation
```

## Plant-Specific QC Considerations

### Inbred Species (rice, soybean, wheat, etc.)
```bash
# Do NOT filter by HWE -- inbred species violate HWE assumption
plink2 --bfile step2_mind \
  --maf 0.01 \
  --make-bed --out cleaned_inbred

# Instead of HWE: filter by heterozygosity
# Remove samples with excessive heterozygosity (potential contamination)
plink2 --bfile cleaned_inbred \
  --het \
  --out het_check

# Then exclude outliers (>3 SD from mean heterozygosity)
```

### Polyploid Species
```bash
# Standard PLINK assumes diploid; be cautious
# Consider using polyploid-aware tools for QC
# UPdog, fitPoly, or polyRAD for initial filtering
# Then convert to pseudo-diploid for standard imputation tools
```

## QC Checkpoints

Record these metrics before and after each filtering step:

| Metric | Command | Expected After QC |
|--------|---------|-------------------|
| Total variants | `wc -l input.bim` | Reported in log |
| Total samples | `wc -l input.fam` | Reported in log |
| Sample retention | After `--mind` | > 90% |
| Variant retention | After all filters | > 50% of input |
| Mean MAF | `plink2 --freq` | Report per chromosome |
| Het rate (inbred) | `plink2 --het` | Flag outliers |

## Strand and Allele Consistency

Before imputation, ensure target and reference alleles match:

```bash
# Check allele coding
plink2 --bfile cleaned \
  --freq --out freq_check

# Flip strands if needed
plink2 --bfile cleaned \
  --flip strand_issue_snps.txt \
  --make-bed --out allele_fixed

# Check against reference
# Use scripts like HRC-1000G-check-bim.pl (adapt for plant references)
perl check_bim_against_reference.pl \
  -b cleaned.bim \
  -r reference_alleles.txt \
  -h
```

## Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| Too few variants after QC | Overly strict thresholds | Relax `--geno` to 0.15-0.20, `--maf` to 0.005 |
| Many samples removed | Poor quality data | Check sample quality; may need re-genotyping |
| HWE filter removes too many | Inbred species | Skip HWE; inbred species violate HWE |
| Chromosome naming mismatch | Different naming conventions | Standardize: `chr1` vs `1`, number-only convention |
| Reference allele mismatch | Strand or coding difference | Run allele checking and alignment script |

## Report Template

After QC, generate a report with:

1. Input: N samples, M variants
2. After `--geno --mind`: N' samples, M' variants
3. After `--maf`: N' samples, M'' variants
4. After `--hwe` (if applicable): N' samples, M''' variants
5. Final: N_final samples, M_final variants
6. Variant retention rate: M_final / M * 100%
7. Any samples flagged for excess heterozygosity or relatedness
