# MetaBAT2

**Goal:** Metagenome binning — group assembled contigs into metagenome-assembled genomes (MAGs)
**Best for:** Complex microbial communities with sufficient sequencing depth (contig coverage >5x)

## Prerequisites

- MetaBAT2 (https://bitbucket.org/berkeleylab/metabat)
- Assembled contigs (FASTA)
- Sorted BAM files for coverage estimation
- `jgi_summarize_bam_contig_depths` for depth calculation

## Basic Usage

```bash
# Step 1: Map reads to assembly for coverage
bowtie2-build assembly.fa assembly_index

for sample in sample*.fq.gz; do
  base=$(basename $sample .fq.gz)
  bowtie2 -x assembly_index \
    -U $sample \
    -S ${base}.sam \
    -p 16
  samtools view -bS ${base}.sam | samtools sort -o ${base}.bam
done

# Step 2: Calculate contig depths
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam

# Step 3: Run MetaBAT2
metabat2 -i assembly.fa \
  -a depth.txt \
  -o bins/bin \
  -t 16 \
  --minContig 1500 \
  --minCV 0.5
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `--minContig` | 1500-2500 | Exclude short plant/host contig fragments |
| `--minCV` | 0.5-0.7 | Minimum coefficient of variation; higher = stricter |
| `--minCVSum` | 1.0 | Minimum total CV across samples |
| `-t` | 16-32 | Thread count |
| `--maxEdges` | 200 | Max edges per node in graph |
| `--pValue` | 0.1 | Higher = more sensitive but more false bins |

## Plant-Specific Notes

- For endosphere samples: use higher `--minContig` (2500) to exclude plant DNA fragments
- For soil samples: use lower `--minContig` (1500) to capture more microbial diversity
- Multiple samples from the same compartment improve differential abundance binning
- If MAG count is <5: likely host removal insufficient, or sequencing depth too low

## Consensus Binning with DASTool

```bash
# Run multiple binners
metabat2 -i assembly.fa -a depth.txt -o metabat_bins/bin
run_concoct.py --coverage_file depth.txt --composition_file assembly.fa -o concoct_bins/
run_maxbin.pl -contig assembly.fa -abund depth.txt -out maxbin_bins/bin

# Integrate with DASTool
DAS_Tool -i metabat_contigs.tsv,concoct_contigs.tsv,maxbin_contigs.tsv \
  -l metabat,concoct,maxbin \
  -c assembly.fa \
  -o dastool_output/ \
  --write_bins 1 \
  --search_engine diamond \
  -t 16
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Too few bins" | Low sequencing depth | Increase depth or lower `--minCV` |
| "Many bins with single contig" | Assembly too fragmented | Check assembly N50; use `--minContig` |
| "Empty bins" | No contigs passed filters | Lower `--minContig` or check depth file format |
| Segment fault (large dataset) | Memory exhausted on large assemblies | Split by contig size or increase RAM |
