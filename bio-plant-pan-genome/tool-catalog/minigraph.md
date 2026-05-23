# Minigraph

**Goal:** Build a sequence skeleton pan-genome graph from multiple genome assemblies
**Best for:** <20 plant genomes, focus on structural variants (SVs) > 50bp

## Prerequisites

- Minigraph (https://github.com/lh3/minigraph)
- FASTA files for each genome assembly
- Enough memory: ~2GB per 100Mb of genome (typical plant: 8-32 GB)

## Basic Usage

```bash
# Build graph from multiple genomes
minigraph -x ggs -t 16 \
  reference.fa \
  genome2.fa genome3.fa genome4.fa \
  > pan_graph.gfa

# Call SVs from the graph
minigraph -x ggs -t 16 --call \
  reference.fa genome2.fa genome3.fa \
  > sv_calls.bed

# Convert GFA to rGFA for downstream tools
minigraph -x ggs -t 16 --rgfa \
  reference.fa genome2.fa \
  > pan_graph.rgfa
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `-x ggs` | ggs (default) | Optimized for genome-to-genome mapping |
| `-t` | 16-32 | Thread count, scale with available CPUs |
| `--call` | optional | Output SV calls in BED format |
| `--rgfa` | optional | Output rGFA format (reference-guided) |
| `-k` | 19 (default) | K-mer size for minimizer; larger = more specific for large genomes |
| `-w` | 10 (default) | Minimizer window size |

## Plant-Specific Notes

- For polyploid genomes, process each subgenome separately
- Minigraph captures SVs >= 50bp but misses SNPs and small indels
- Large plant genomes (wheat ~17Gb): use `-k 29 -w 15` for higher specificity
- Output GFA can be visualized with Bandage or ODGI

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Out of memory | Large genomes | Split by chromosome or reduce threads |
| "No alignment found" | Too divergent genomes | Lower `-k` value to 15 |
| Empty GFA | Input FASTA issue | Check FASTA headers — no spaces allowed |
| Segmentation fault | Corrupt FASTA | Validate with `seqtk seq -A input.fa > /dev/null` |
