# Beagle 5 -- Genotype Imputation and Phasing

**Goal:** Impute missing genotypes and phase haplotypes without requiring an external reference panel
**Best for:** Small-to-medium plant datasets (n<500, SNP<50K), polyploids, non-model species

## Prerequisites

- Java 8+ (Beagle is a Java JAR)
- Input: VCF (gzipped or not) with GT field
- Optional: reference VCF (can use internal reference if omitted)
- Download from: https://faculty.washington.edu/browning/beagle/beagle.html

## Basic Usage

### Without Reference Panel (internal imputation)

```bash
java -Xmx8g -jar beagle.22Jul22.46e.jar \
  gt=input.vcf.gz \
  out=output_prefix \
  ne=1000 \
  window=40.0 \
  nthreads=4
```

### With Reference Panel

```bash
java -Xmx16g -jar beagle.22Jul22.46e.jar \
  gt=input.vcf.gz \
  ref=reference_panel.vcf.gz \
  out=output_prefix \
  impute=true \
  ne=1000 \
  nthreads=8
```

### Polyploid Mode (experimental)

```bash
java -Xmx16g -jar beagle.22Jul22.46e.jar \
  gt=tetraploid_input.vcf.gz \
  out=tetraploid_output \
  ne=500 \
  nthreads=4
```

## Key Parameters

| Parameter | Default | Plant Recommendation |
|-----------|---------|---------------------|
| `ne` | 1000000 | 100-1000 for breeding populations; 200-500 for inbred species |
| `window` | 40.0 cM | 40-60 for inbred species (longer LD); 20-40 for outcross species |
| `nthreads` | 1 | Match available CPU cores |
| `impute` | true | Keep true for imputation; false for phasing only |
| `overlap` | 3.0 cM | Default fine; increase to 5.0 for low-density markers |
| `burnin-its` | 5 | 5-10; increase for higher accuracy at cost of runtime |
| `phase-its` | 5 | 5-10; same trade-off as burnin |

## Output Files

- `output_prefix.vcf.gz` -- Imputed and phased VCF
- `output_prefix.log` -- Run log with timing and memory usage
- DR2 (dosage R-squared) is embedded in the VCF INFO field as `DR2=`

## Plant Relevance

- **Best tool for non-model species**: No reference panel needed; works with any species
- **Polyploid support**: Experimental but better than IMPUTE2/Minimac4 which assume diploid
- **Inbred species tip**: Lower `ne` (200-500) and wider `window` (50-60 cM) to account for extended LD
- **Outcross species tip**: Keep `ne` higher (500-2000) and narrower `window` (20-40 cM)

## Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| `OutOfMemoryError` | Java heap too small | Increase `-Xmx` (e.g., `-Xmx32g`) |
| `No GT field` | VCF lacks genotype data | Check input VCF; add GT with `bcftools` |
| `Too few markers` | < 100 variants in region | Merge with nearby regions or skip |
| Very slow runtime | Too many samples/markers | Reduce `ne`, `burnin-its`, `phase-its`; split by chromosome |
| DR2 < 0.3 for most variants | Insufficient LD information | Increase sample size or marker density; consider using reference panel |

## Runtime Estimates

| Data Size | Typical Runtime |
|-----------|----------------|
| 1K SNPs, 100 samples | < 1 minute |
| 10K SNPs, 300 samples | 5-15 minutes |
| 50K SNPs, 500 samples | 30-90 minutes |
| 500K SNPs, 1000 samples | 4-12 hours |

## Validation

To assess imputation accuracy, mask 5-10% of known genotypes and compare:

```bash
# Mask random genotypes
bcftools +setGT input.vcf.gz -- -t q -n . -r mask.bed > masked.vcf

# Impute masked file
java -jar beagle.22Jul22.46e.jar gt=masked.vcf out=imputed

# Compare (requires custom script or bcftools plugin)
bcftools +fill-tags imputed.vcf.gz -- -t all | \
  bcftools query -f '%CHROM\t%POS\t%DR2\n' > accuracy_report.txt
```
