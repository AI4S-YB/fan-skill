# BSA Analysis Primer

## Analysis Overview

Bulk Segregant Analysis (BSA) is a rapid gene mapping technique that identifies genomic regions associated with a trait by comparing allele frequencies between two pooled DNA samples with contrasting phenotypes.

## Key Methods

### 1. MutMap
- **Design**: Mutant × Wild-type cross, F2 mutant pool
- **Signal**: SNP-index approaching 1.0 for causal mutations
- **Best for**: Identifying causal mutations in mutants

### 2. QTL-seq
- **Design**: Two pools from biparental cross (extreme phenotypes)
- **Signal**: ΔSNP-index deviation from 0
- **Best for**: QTL mapping in biparental populations

### 3. BSR-seq
- **Design**: RNA-seq based BSA
- **Signal**: Combines SNP-index with expression differences
- **Best for**: When DNA-seq is not available; expression + mapping combined

## Input Requirements

| Input | Description | Required |
|-------|-------------|----------|
| Bulk pools FASTQ | Sequencing data for 2+ pools | Yes |
| Reference genome | FASTA format | Yes |
| GFF annotation | Gene annotations | Recommended |
| Parent data | Parental sequencing | For MutMap |
| Population info | F2/RIL/BC generation | Yes |

## Critical Parameters

### Pool Size
| Population | Minimum | Recommended |
|------------|---------|-------------|
| F2 | 20 | 50-100 |
| RIL | 20 | 30-50 |
| BC | 30 | 50-100 |

### Sequencing Depth
| Depth | Quality |
|-------|---------|
| < 10x | Insufficient |
| 10-20x | Marginal |
| 20-50x | Good |
| > 50x | Excellent |

### Sliding Window
| Parameter | Typical Value |
|-----------|---------------|
| Window size | 1 Mb |
| Step size | 100 kb |
| Min SNPs/window | 10 |

## Expected Outputs

| Output | Description |
|--------|-------------|
| SNP-index profile | Per-SNP allele frequency |
| Sliding window plot | Smoothed signal |
| Confidence intervals | Statistical significance |
| Candidate intervals | Genomic regions |
| Candidate genes | Annotated genes in regions |

## Common Pitfalls

1. **Insufficient pool size** → False negatives
2. **Low sequencing depth** → High variance in SNP-index
3. **Population structure** → Spurious signals
4. **Recombination rate variation** → Signal spread varies
5. **Multiple linked QTL** → Complex signal patterns

## Tools Overview

| Tool | Purpose |
|------|---------|
| BWA-MEM | Read alignment |
| GATK HaplotypeCaller | Variant calling |
| bcftools | Variant processing |
| QTL-seq pipeline | End-to-end analysis |
| MutMap pipeline | MutMap specific |
| SnpEff | Variant annotation |
