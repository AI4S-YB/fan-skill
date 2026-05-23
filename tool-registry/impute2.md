# IMPUTE2 -- Reference-Panel-Based Genotype Imputation

**Goal:** High-accuracy imputation using a haplotype reference panel
**Best for:** Species with well-characterized reference panels (e.g., rice 3K, maize HapMap)

## Prerequisites

- IMPUTE2 software: https://mathgen.stats.ox.ac.uk/impute/impute_v2.html
- SHAPEIT (for pre-phasing): https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html
- Haplotype reference panel in IMPUTE2 format (.hap / .legend / .map)
- Genetic map (recombination rate per marker, cM)
- Input: phased or unphased VCF/PLINK

## Workflow

### Step 1: Pre-phasing with SHAPEIT

```bash
shapeit \
  --input-vcf input.vcf.gz \
  --input-map genetic_map.txt \
  --output-max phased \
  --thread 4
```

### Step 2: Imputation with IMPUTE2

```bash
impute2 \
  -m genetic_map.txt \
  -h reference.hap.gz \
  -l reference.legend.gz \
  -g phased.gen.gz \
  -Ne 500 \
  -buffer 250 \
  -burnin 10 \
  -iter 30 \
  -k_hap 500 \
  -o output_prefix \
  -int 1 50000000
```

### Step 3: Convert output to VCF

```bash
# IMPUTE2 outputs .gen format; convert with GTOOL or bcftools
gtool -G --g output_prefix.gen --s sample_file.sample --ped output.ped --map output.map
plink --file output --recode vcf --out output
```

## Key Parameters

| Parameter | Default | Plant Recommendation |
|-----------|---------|---------------------|
| `-Ne` | 20000 | 100-1000 for plant breeding populations; 200-500 for inbreds |
| `-k_hap` | 80 | 500-1000 for better accuracy (costs runtime) |
| `-buffer` | 250 kb | 250-500 kb for species with long LD (inbred crops) |
| `-burnin` | 10 | 10-20 for stable convergence |
| `-iter` | 30 | 30-50; more iterations improve accuracy at cost of runtime |
| `-int` | whole chr | Split chromosomes into 5-10 Mb chunks for parallel execution |

## Output Files

- `output_prefix` -- Imputed genotypes (.gen format with dosage probabilities)
- `output_prefix_info` -- Per-variant info score (0-1, higher is better)
- `output_prefix_summary` -- Run summary

## Plant Relevance

- **Requires reference panel**: Only suitable when a species-specific reference exists
- **Known plant references**: rice 3K, maize HapMap3, soybean 302 accessions, wheat 10+ genomes
- **Inbred species**: Can use smaller reference panels due to extended LD
- **Phasing requirement**: SHAPEIT pre-phasing adds an extra step and tool dependency

## Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| `Segmentation fault` | Memory limit exceeded | Split by chromosome; use smaller `-k_hap` |
| `Reference mismatch` | Target/ref allele coding inconsistent | Standardize alleles (flip/swap) before imputation |
| `Low info score (<0.4)` | Reference panel poorly matched | Check population match; consider Beagle instead |
| `SHAPEIT fails on polyploid` | SHAPEIT is diploid-only | Use Beagle for polyploid species |
| Very slow per chunk | `-k_hap` too large | Reduce to 300-500; test runtime on small region first |

## Runtime Estimates

| Data Size | Typical Runtime (per chromosome) |
|-----------|----------------------------------|
| 10K SNPs, 200 samples, 500 ref haps | 1-3 hours |
| 50K SNPs, 500 samples, 1000 ref haps | 6-24 hours |
| 500K SNPs, 2000 samples | Several days; use cluster parallelization |

## Info Score Interpretation

- **info > 0.9**: Excellent -- suitable for GWAS and fine mapping
- **info 0.7-0.9**: Good -- acceptable for most analyses
- **info 0.4-0.7**: Marginal -- use with caution; avoid for low-MAF variants
- **info < 0.4**: Poor -- exclude from downstream analysis
- **info < 0.3**: Unreliable -- always exclude
