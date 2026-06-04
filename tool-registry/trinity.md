# Trinity: De Novo Transcriptome Assembly

## Tool Overview

**Tool ID**: `trinity`
**Category**: Transcriptome Assembly
**Purpose**: De novo assembly of RNA-seq data without reference genome
**Web**: https://github.com/trinityrnaseq/trinityrnaseq
**Version**: Trinity v2.15+

## When to Use

- Non-model plant species without reference genome
- Transcript discovery in unannotated species
- Alternative splicing analysis without genome
- lncRNA discovery

## Key Parameters

| Parameter | Default | Recommended for Plants | Description |
|-----------|---------|------------------------|-------------|
| `--min_contig_length` | 200 | **200** | Minimum assembled contig length |
| `--max_memory` | - | **100G-200G** | Maximum RAM for Jellyfish |
| `--CPU` | 1 | **16-32** | Number of threads |
| `--path_reinforcement_distance` | 75 | **75-150** | Butterfly path reinforcement |
| `--min_kmer_cov` | 1 | **2** | Minimum k-mer coverage |
| `--full_cleanup` | false | true | Clean intermediate files |

### Critical Parameters for Plants

#### 1. Memory Requirements

Plant transcriptomes vary widely:
- Small genome (<1Gb): 50-100G RAM
- Medium genome (1-2Gb): 100-200G RAM
- Large genome (>2Gb): 200-500G RAM

```bash
# Always specify max_memory for plants
Trinity --seqType fq \
  --left R1.fq.gz --right R2.fq.gz \
  --max_memory 200G \
  --CPU 32 \
  --output trinity_out
```

#### 2. Path Reinforcement Distance

For plants with complex gene families:

```bash
# Increase for large gene families
--path_reinforcement_distance 150
```

## Command Examples

### Standard Assembly

```bash
# Basic Trinity assembly
Trinity --seqType fq \
  --left sample_R1.fq.gz \
  --right sample_R2.fq.gz \
  --max_memory 100G \
  --CPU 16 \
  --output trinity_assembly
```

### Multiple Samples (Combined Assembly)

```bash
# Create samples file
# samples.txt:
# cond_A	rep1	/path/to/A_rep1_R1.fq.gz	/path/to/A_rep1_R2.fq.gz
# cond_A	rep2	/path/to/A_rep2_R1.fq.gz	/path/to/A_rep2_R2.fq.gz
# cond_B	rep1	/path/to/B_rep1_R1.fq.gz	/path/to/B_rep1_R2.fq.gz

Trinity --seqType fq \
  --samples_file samples.txt \
  --max_memory 200G \
  --CPU 32 \
  --output trinity_combined
```

### With Normalization (Large Datasets)

```bash
# Normalize reads to reduce memory
Trinity --seqType fq \
  --left R1.fq.gz --right R2.fq.gz \
  --normalize_by_read_set \
  --max_memory 100G \
  --CPU 16 \
  --output trinity_normalized
```

### Strand-Specific Libraries

```bash
# For strand-specific RNA-seq
Trinity --seqType fq \
  --left R1.fq.gz --right R2.fq.gz \
  --SS_lib_type RF \
  --max_memory 100G \
  --CPU 16 \
  --output trinity_stranded
```

**Strand types**:
- `RF`: First-strand (dUTP method) - **most common for plants**
- `FR`: Second-strand
- `F`: Unpaired, forward
- `R`: Unpaired, reverse

## Output Files

| File | Description |
|------|-------------|
| `Trinity.fasta` | Final assembled transcripts |
| `Trinity.fasta.gene_trans_map` | Gene-transcript mapping |
| `Trinity.timing` | Runtime statistics |
| `TrinityStats.txt` | Assembly statistics |

## Assembly Quality Assessment

### N50 and Statistics

```bash
# Use TrinityStats.pl
$TRINITY_HOME/util/TrinityStats.pl trinity_out/Trinity.fasta

# Output example:
# Total trinity 'genes':  50000
# Total trinity transcripts:  75000
# Percent GC: 45.2
# Contig N50: 1500
# Median contig length: 800
```

### BUSCO Completeness

```bash
# Assess completeness with plant BUSCOs
busco -i Trinity.fasta \
  -l embryophyta_odb10 \
  -m transcriptome \
  -o busco_output \
  -c 16
```

**Expected for good assembly**:
- Complete BUSCOs: > 80%
- Fragmented: < 10%
- Missing: < 10%

### TransRate

```bash
# Comprehensive quality assessment
transrate --assembly Trinity.fasta \
  --left R1.fq.gz --right R2.fq.gz \
  --threads 16 \
  --output transrate_out
```

## Plant-Specific Considerations

### 1. High Heterozygosity

For highly heterozygous plants:

```bash
# May produce redundant haplotypes
# Solution: Post-assembly deduplication
cd-hit-est -i Trinity.fasta -o Trinity_nr.fasta -c 0.95
```

### 2. Polyploid Species

Polyploids present challenges:
- Homeologs may assemble as separate transcripts
- Higher transcript counts expected
- Consider subgenome-specific assembly

### 3. Gene Family Expansion

Plants have expanded gene families:

```bash
# Expect higher transcript counts
# 50,000-150,000 transcripts common for plants
# Check for gene family representation
```

## Integration with Downstream Analysis

### Expression Quantification (RSEM)

```bash
# Build RSEM reference
$TRINITY_HOME/util/RSEM/reference/rsem-prepare-reference \
  Trinity.fasta trinity_rsem

# Quantify expression
rsem-calculate-expression --paired-end \
  --num-threads 8 \
  sample_R1.fq.gz sample_R2.fq.gz \
  trinity_rsem sample_name
```

### Expression Quantification (Salmon)

```bash
# Build Salmon index
salmon index -t Trinity.fasta -i trinity_index --type quasi -k 31

# Quantify
salmon quant -i trinity_index \
  -l A \
  -1 sample_R1.fq.gz \
  -2 sample_R2.fq.gz \
  -p 8 \
  -o sample_quant
```

### Functional Annotation (Trinotate)

```bash
# Full annotation pipeline
# 1. Predict ORFs
TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta

# 2. Run BLAST and hmmscan
blastp -query longest_orfs.pep -db swissprot -out blast.out
hmmscan --cpu 8 Pfam-A.hmm longest_orfs.pep > pfam.out

# 3. Generate annotation
$TRINITY_HOME/trinotate/Trinotate Trinity.fasta \
  --gene_trans_map Trinity.fasta.gene_trans_map \
  --transdecoder_pep longest_orfs.pep \
  --pfam_domain pfam.out \
  --blastp blast.out \
  > trinotate_annotation.tsv
```

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| Out of memory | Insufficient RAM | Increase `--max_memory`, use normalization |
| Too many transcripts | Heterozygosity, contamination | Use CD-HIT to deduplicate |
| Low N50 | Low coverage, poor quality | Check sequencing depth (>30M reads) |
| Missing BUSCOs | Incomplete assembly | Increase coverage, check RNA quality |
| Long runtime | Large dataset | Increase `--CPU`, use grid computing |

## References

- Grabherr et al. (2011) Trinity: reconstructing a full-length transcriptome without a genome
- GitHub: https://github.com/trinityrnaseq/trinityrnaseq

## Related Tools

- **TransDecoder**: ORF prediction
- **Trinotate**: Functional annotation
- **RSEM**: Expression quantification
- **Salmon**: Fast quantification
- **BUSCO**: Completeness assessment
- **CD-HIT**: Sequence clustering
