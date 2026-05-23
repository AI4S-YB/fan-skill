# Minimac4 -- Efficient Genotype Imputation for Large Samples

**Goal:** Fast genotype imputation using state-space reduction
**Best for:** Large sample sizes (>=500) with M3VCF reference panels

## Prerequisites

- Minimac4: https://github.com/statgen/Minimac4
- Reference panel in M3VCF format (convert from VCF using Minimac3)
- Input: VCF (phased or unphased)

## Workflow

### Step 1: Prepare M3VCF reference

```bash
# Convert VCF reference to M3VCF using Minimac3
Minimac3 \
  --refHaps reference.vcf.gz \
  --processReference \
  --prefix reference \
  --cpus 8
```

### Step 2: Run Minimac4 imputation

```bash
minimac4 \
  --refHaps reference.m3vcf.gz \
  --haps target.vcf.gz \
  --prefix output_prefix \
  --cpus 8 \
  --minRatio 0.01 \
  --meta
```

### Alternative: One-step with VCF reference (Minimac4 >= 1.0)

```bash
minimac4 \
  --refHaps reference.vcf.gz \
  --haps target.vcf.gz \
  --prefix output_prefix \
  --cpus 8
```

## Key Parameters

| Parameter | Default | Plant Recommendation |
|-----------|---------|---------------------|
| `--cpus` | 1 | Match available cores; Minimac4 scales well |
| `--minRatio` | 0.01 | 0.001 for more variants; 0.01 for common variants only |
| `--meta` | off | Enable for detailed per-variant metrics |
| `--probThreshold` | 0.01 | 0.001 for more precise genotype probabilities |
| `--diffThreshold` | 0.01 | Adjust for convergence in rare variants |
| `--topThreshold` | - | Limit number of reference haplotypes considered |
| `--chunkLengthMb` | 20 | 20-50 Mb; smaller chunks for better memory usage |

## Output Files

- `output_prefix.dose.vcf.gz` -- Dosage VCF (0-2 continuous values)
- `output_prefix.info` -- Per-variant imputation quality (R², MAF, etc.)
- `output_prefix.empiricalDose.vcf.gz` -- Empirical dosage (with `--meta`)

## Plant Relevance

- **Best for large breeding datasets**: Genomic selection panels with 500+ samples
- **Requires M3VCF reference**: Extra step to prepare reference, but amortized over many imputations
- **Good for species with established reference panels**: maize HapMap, rice 3K, soybean
- **Speed advantage**: Imputes 1000+ samples in hours, not days

## Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| `M3VCF index missing` | Reference not indexed | Re-run with `--processReference` |
| `Allele mismatch` | Target/ref alleles differ | Pre-process target VCF to match reference alleles |
| `Out of memory` | Chunk too large | Reduce `--chunkLengthMb` to 10-15 |
| `Empty output` | No overlapping variants | Check chromosome naming consistency |
| `Low R² for many variants` | Reference panel too distant | Consider Beagle internal imputation instead |

## Info Score (Empirical R²) Interpretation

Minimac4 reports R² (estimated squared correlation between imputed and true genotypes):

- **R² > 0.8**: High confidence -- suitable for all downstream analyses
- **R² 0.5-0.8**: Moderate -- usable for most purposes, flag low-MAF variants
- **R² 0.3-0.5**: Low -- use with caution; exclude for GWAS
- **R² < 0.3**: Very low -- exclude from analysis

## Runtime Estimates

| Data Size | Typical Runtime |
|-----------|----------------|
| 50K SNPs, 500 samples, 1M ref | 1-3 hours |
| 500K SNPs, 1000 samples, 10M ref | 6-12 hours |
| 1M+ SNPs, 5000+ samples | 1-3 days (with parallelization per chromosome) |

## Minimac4 vs. IMPUTE2

| Aspect | Minimac4 | IMPUTE2 |
|--------|----------|---------|
| Speed | 10-50x faster | Slower |
| Memory | Lower | Higher |
| Reference format | M3VCF | .hap/.legend |
| Pre-phasing | Optional (better with) | Required |
| Plant panels | Limited M3VCF | More .hap/.legend |
| Output | Dosage VCF | .gen (needs conversion) |
| Community | Growing | Established but aging |
