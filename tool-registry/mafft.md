# MAFFT: Multiple Sequence Alignment

## Tool Overview

**Tool ID**: `mafft`
**Category**: Multiple Sequence Alignment
**Purpose**: Multiple sequence alignment for proteins and nucleotides
**Source**: https://mafft.cbrc.jp/alignment/software/
**Version**: MAFFT 7+ (recommended)

## When to Use

- Gene family sequence alignment
- Phylogenetic tree input preparation
- HMM profile building input
- Codon alignment for Ka/Ks analysis

## Alignment Strategies

| Strategy | Command | Accuracy | Speed | Best For |
|----------|---------|----------|-------|----------|
| **L-INS-i** | `--localpair --maxiterate 1000` | Highest | Slow | <50 sequences, publication quality |
| **G-INS-i** | `--globalpair --maxiterate 1000` | High | Medium | Global alignment, similar length |
| **E-INS-i** | `--genafpair --maxiterate 1000` | High | Slow | Sequences with large gaps |
| **FFT-NS-2** | `--retree 2` | Good | Fast | ≥50 sequences |
| **FFT-NS-1** | `--retree 1` | Fair | Fastest | >500 sequences, quick check |
| **Auto** | `--auto` | Adaptive | Variable | Let MAFFT decide |

## Key Parameters

| Parameter | Default | Recommended | Description |
|-----------|---------|-------------|-------------|
| `--auto` | - | **Use often** | Automatic strategy selection |
| `--localpair` | - | L-INS-i mode | Local pairwise alignment |
| `--globalpair` | - | G-INS-i mode | Global pairwise alignment |
| `--maxiterate` | 0 | **1000** | Iterative refinement cycles |
| `--retree` | 2 | 1-2 | Guide tree iterations |
| `--thread` | 1 | **8** | Number of threads |
| `--adjustdirection` | - | Useful | Adjust strand orientation |
| `--add` | - | - | Add sequences to existing alignment |

## Command Examples

### Small Alignment (<50 sequences, High Quality)

```bash
# L-INS-i: Most accurate for small sets
mafft --localpair --maxiterate 1000 --thread 8 \
  input_sequences.fasta > aligned.fasta

# Short form
mafft --linsi --thread 8 input.fasta > aligned.fasta
```

### Medium Alignment (50-200 sequences)

```bash
# FFT-NS-2: Balance of speed and accuracy
mafft --retree 2 --maxiterate 2 --thread 8 \
  input_sequences.fasta > aligned.fasta

# Or use auto
mafft --auto --thread 8 input.fasta > aligned.fasta
```

### Large Alignment (>200 sequences)

```bash
# FFT-NS-1: Fast for large sets
mafft --retree 1 --thread 16 \
  input_sequences.fasta > aligned.fasta

# Or use auto (will select fast method)
mafft --auto --thread 16 input.fasta > aligned.fasta
```

### Very Large Alignment (>1000 sequences)

```bash
# Use parttree algorithm for huge alignments
mafft --parttree --thread 16 input.fasta > aligned.fasta
```

### Adding Sequences to Existing Alignment

```bash
# Add new sequences to pre-aligned dataset
mafft --add new_sequences.fasta \
  --thread 8 \
  existing_alignment.fasta > combined_alignment.fasta
```

### Codon Alignment for Ka/Ks

```bash
# For Ka/Ks analysis, align proteins first, then map to CDS

# Step 1: Align protein sequences
mafft --auto --thread 8 proteins.fasta > aligned_proteins.fasta

# Step 2: Use PAL2NAL for codon alignment
pal2nal.pl aligned_proteins.fasta cds_sequences.fasta \
  -output fasta > codon_aligned.fasta
```

## Plant-Specific Considerations

### 1. Multi-Domain Proteins

Plant proteins often have multiple domains with varying conservation:

```bash
# Use E-INS-i for proteins with large insertions/deletions
mafft --genafpair --maxiterate 1000 --thread 8 \
  multi_domain_proteins.fasta > aligned.fasta
```

### 2. Polyploid Gene Families

For gene families with many similar homeologs:

```bash
# High accuracy for distinguishing similar sequences
mafft --localpair --maxiterate 1000 --thread 8 \
  --adjustdirection \
  homeologs.fasta > aligned.fasta
```

### 3. Large Gene Families

Plant gene families can have 100+ members:

```bash
# For families >100 members
mafft --auto --thread 16 large_family.fasta > aligned.fasta

# If quality is critical, use FFT-NS-2 with more iterations
mafft --retree 2 --maxiterate 100 --thread 16 \
  large_family.fasta > aligned.fasta
```

## Integration with Gene Family Workflow

### Phylogenetic Tree Input

```bash
# Complete workflow: Alignment → Trimming → Tree

# 1. Align sequences
mafft --auto --thread 8 family_proteins.fasta > aligned.fasta

# 2. Trim poorly aligned regions
trimal -in aligned.fasta -out aligned_trimmed.fasta -automated1

# 3. Build phylogenetic tree
iqtree2 -s aligned_trimmed.fasta -m MFP -B 1000 --alrt 1000 -T 8
```

### HMM Profile Building

```bash
# Build HMM profile from alignment

# 1. Align known family members
mafft --linsi --thread 8 known_members.fasta > aligned.fasta

# 2. Build HMM
hmmbuild family.hmm aligned.fasta

# 3. Optionally calibrate
hmmcalibrate family.hmm
```

### Combined Pipeline

```bash
#!/bin/bash
# gene_family_alignment.sh

INPUT=$1
OUTPUT_DIR=$2
THREADS=${3:-8}

mkdir -p $OUTPUT_DIR

# Detect sequence count
SEQ_COUNT=$(grep -c "^>" $INPUT)

# Choose strategy based on count
if [ $SEQ_COUNT -lt 50 ]; then
  echo "Using L-INS-i (accurate) for $SEQ_COUNT sequences"
  mafft --localpair --maxiterate 1000 --thread $THREADS \
    $INPUT > $OUTPUT_DIR/aligned.fasta
elif [ $SEQ_COUNT -lt 200 ]; then
  echo "Using FFT-NS-2 (balanced) for $SEQ_COUNT sequences"
  mafft --retree 2 --maxiterate 2 --thread $THREADS \
    $INPUT > $OUTPUT_DIR/aligned.fasta
else
  echo "Using FFT-NS-1 (fast) for $SEQ_COUNT sequences"
  mafft --retree 1 --thread $THREADS \
    $INPUT > $OUTPUT_DIR/aligned.fasta
fi

# Trim alignment
trimal -in $OUTPUT_DIR/aligned.fasta \
  -out $OUTPUT_DIR/aligned_trimmed.fasta \
  -automated1

echo "Alignment complete: $OUTPUT_DIR/aligned_trimmed.fasta"
```

## Alignment Quality Assessment

```bash
# Check alignment quality with ALISCORE
aliscore $OUTPUT_DIR/aligned_trimmed.fasta

# Or use trimAl's check
trimal -in aligned.fasta -gt 0.5 -cons 60 -out filtered.fasta
```

## Output Format

MAFFT outputs standard FASTA alignment:

```
>Gene001
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGC---
>Gene002
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG
>Gene003
ATGCGTACGTAGCTAGCTAGCTAGCTAG------C
```

## Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| Slow for large sets | Wrong strategy | Use `--auto` or FFT-NS methods |
| Poor alignment | Divergent sequences | Try L-INS-i or E-INS-i |
| Memory error | Too many sequences | Use `--parttree` |
| Misaligned region | Repeat domains | Use E-INS-i for large gaps |
| Strand direction wrong | Mixed orientation | Use `--adjustdirection` |

## Performance Comparison

| Sequences | Method | Time (approx) | Memory |
|-----------|--------|---------------|--------|
| 10 | L-INS-i | seconds | low |
| 50 | L-INS-i | minutes | medium |
| 50 | FFT-NS-2 | seconds | low |
| 200 | FFT-NS-2 | minutes | medium |
| 1000 | FFT-NS-1 | minutes | medium |
| 10000 | Parttree | minutes | high |

## References

- Katoh & Standley (2013) MAFFT multiple sequence alignment software
- Katoh et al. (2002) MAFFT: a novel method for rapid multiple sequence alignment
- MAFFT Documentation: https://mafft.cbrc.jp/alignment/software/manual/manual.html

## Related Tools

- **trimAl**: Alignment trimming and quality filtering
- **PAL2NAL**: Codon alignment from protein alignment
- **IQ-TREE**: Phylogenetic tree construction
- **HMMER**: HMM profile building from alignment
- **Clustal Omega**: Alternative aligner
- **MUSCLE**: Alternative aligner (faster for small sets)
