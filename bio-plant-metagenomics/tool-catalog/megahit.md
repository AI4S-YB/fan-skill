# MEGAHIT

**Goal:** Metagenome assembly from short reads using succinct de Bruijn graphs
**Best for:** Medium to large metagenomic datasets (>5 Gbp), complex microbial communities

## Prerequisites

- MEGAHIT (https://github.com/voutcn/megahit)
- Paired-end FASTQ files (quality-trimmed and host-filtered)
- Memory: ~1 byte per k-mer for soil metagenomes (typically 50-200 GB)

## Basic Usage

```bash
# Standard assembly
megahit -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
  -o megahit_output/ \
  -t 32 \
  --k-list 21,29,39,59,79,99,119 \
  --min-contig-len 500

# With memory limit
megahit -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
  -o megahit_output/ \
  -t 32 \
  -m 0.8 \
  --k-list 27,47,67,87,107 \
  --min-contig-len 1000

# Co-assembly (multiple samples)
megahit -1 sampleA_R1.fq,sampleB_R1.fq,sampleC_R1.fq \
  -2 sampleA_R2.fq,sampleB_R2.fq,sampleC_R2.fq \
  -o coassembly_output/ \
  -t 64 \
  --k-min 27 --k-max 127 --k-step 20
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `--k-list` | 21,29,39,59,79,99,119 | Range captures low to high coverage genomes |
| `--min-contig-len` | 500-1000 | Shorter for sensitivity, longer for quality |
| `-m` | 0.8 | Fraction of system memory to use |
| `-t` | 32-64 | Thread count |
| `--k-min/--k-max/--k-step` | 27/127/20 | Alternative to explicit k-list |

## Plant-Specific Notes

- Soil/rhizosphere samples have the highest complexity; use wide k-mer range
- After host removal, microbial read depth drops — use `--min-count 2` for low-depth samples
- Co-assembly across samples from the same compartment increases assembly continuity
- Expect N50 of 1-5 kbp for soil, 5-20 kbp for endosphere (simpler communities)

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Out of memory" | Too many k-mers or reads | Reduce `-m`, narrow k-mer range, or subsample reads |
| "Too few contigs" | Insufficient microbial reads | Check host removal efficiency |
| "All reads filtered" | Quality too strict | Check input FASTQ quality with FastQC |
| Slow assembly | Very large dataset (>50 Gbp) | Use `--no-mercy` for highly complex communities |
