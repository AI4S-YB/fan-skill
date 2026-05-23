# PGGB (PanGenome Graph Builder)

**Goal:** Build a full-variation pan-genome graph via all-vs-all alignment
**Best for:** 20-100 plant genomes where complete variation (SNP + indel + SV) is needed

## Prerequisites

- PGGB (https://github.com/pangenome/pggb)
- Docker/Singularity recommended for dependency management
- FASTA files with standardized headers
- Large memory: 128-512 GB for plant genomes

## Basic Usage

```bash
# Basic PGGB run
pggb -i genome_list.txt \
  -o pggb_output/ \
  -t 32 \
  -p 95 \
  -s 100000 \
  -n 20 \
  -G 20000,200000 \
  -k 19

# With segment/POA length parameters
pggb -i genome_list.txt \
  -o pggb_output/ \
  -t 32 \
  -p 95 \
  -s 50000 \
  -n 10 \
  -K 79

# Chromosome-by-chromosome (for large genomes)
for chr in Chr01 Chr02 Chr03; do
  pggb -i chr_${chr}_list.txt \
    -o pggb_${chr}/ \
    -t 16 \
    -p 95 \
    -s 100000
done
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `-p` | 90-98 | Percent identity for mapping — lower for divergent species |
| `-s` | 50000-200000 | Segment length for pairwise alignment |
| `-n` | 10-20 | Number of haplotypes/genomes |
| `-G` | 20000,200000 | Min/max POA length |
| `-k` | 19-79 | K-mer size; larger for large genomes |
| `-t` | 16-64 | Thread count |

## Plant-Specific Notes

- For polyploids: process each subgenome separately, then merge graphs
- Large genomes (>1Gb): split by chromosome, run PGGB per chromosome
- Divergent plant species: lower `-p` to 85-90
- Rice (400Mb): ~30 min per chromosome for 20 genomes
- Wheat (17Gb): ~6-12 hours per chromosome for 10 genomes

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "wfmash segmentation fault" | Memory exhausted | Reduce `-s`, split chromosomes, or use more RAM |
| "seqwish OOM" | Too many genomes in one run | Reduce genomes per batch or use `-P` for pruning |
| "no alignment produced" | Genomes too divergent | Lower `-p` or pre-align with different parameters |
| "too many nodes" | High diversity | Increase `-G` min POA length or filter low-complexity regions |
