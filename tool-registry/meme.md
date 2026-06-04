# MEME: Motif Discovery and Analysis

## Tool Overview

**Tool ID**: `meme`
**Category**: Motif Analysis
**Purpose**: De novo motif discovery in DNA or protein sequences
**Web**: https://meme-suite.org/
**Version**: MEME Suite 5.x

## When to Use

- Gene family motif analysis
- Transcription factor binding site discovery
- Protein domain identification
- Regulatory element analysis

## Key Parameters

| Parameter | Default | Recommended for Plants | Description |
|-----------|---------|------------------------|-------------|
| `-mod` | zoops | **anr** (plants) | Motif distribution model |
| `-nmotifs` | 1 | **10** | Number of motifs to find |
| `-minw` | 8 | **6** | Minimum motif width |
| `-maxw` | 50 | **50** | Maximum motif width |
| `-objfun` | classic | classic | Objective function |
| `-revcomp` | no | **yes** (DNA) | Search both strands |

### Critical Parameter: `-mod` (Distribution Model)

**Why `-mod anr` for plants?**

| Model | Description | Best For |
|-------|-------------|----------|
| **oops** | One Occurrence Per Sequence | Ubiquitous TF motifs |
| **zoops** | Zero or One Occurrence Per Sequence | Animal TFs, clean data |
| **anr** | Any Number of Repetitions | **Plant proteins**, TFs with repeats |

**Plant-specific reasoning**:
- Plant protein domains often repeat (e.g., WRKY, NBS-LRR)
- Same motif may appear multiple times in one sequence
- `anr` model captures these repeats without penalty

### Recommended Parameter Sets

#### Protein Sequences (Gene Family)

```bash
meme proteins.fasta \
  -mod anr \
  -nmotifs 10 \
  -minw 6 \
  -maxw 50 \
  -protein \
  -oc meme_output \
  -cp 8
```

#### DNA Sequences (Promoters)

```bash
meme promoters.fasta \
  -mod anr \
  -nmotifs 15 \
  -minw 6 \
  -maxw 20 \
  -dna \
  -revcomp \
  -oc meme_output \
  -cp 8
```

## Output Files

| File | Content |
|------|---------|
| `meme.html` | Interactive results page |
| `meme.txt` | Text format results |
| `meme.xml` | XML format for downstream tools |
| `meme.erect.txt` | Erect format for visualization |

## Motif Discovery Workflow

### Step 1: Prepare Input Sequences

```bash
# For protein sequences - use full-length or specific domains
# FASTA format required

# Remove sequences with >50% gaps if using alignment
seqkit seq -g aligned_proteins.fasta > clean_proteins.fasta
```

### Step 2: Run MEME

```bash
# Standard gene family analysis
meme gene_family_proteins.fasta \
  -mod anr \
  -nmotifs 10 \
  -minw 6 \
  -maxw 50 \
  -protein \
  -oc meme_results \
  -cp 8
```

### Step 3: Motif Annotation with TOMTOM

```bash
# Compare discovered motifs to known databases
tomtom -oc tomtom_results \
  -db JASPAR_CORE_2022.meme \
  -db $MEME_DB/motif_databases/PLANTTFDB.meme \
  meme_results/meme.txt
```

### Step 4: Motif Visualization

```bash
# Generate sequence logo
# Use online tool: http://weblogo.threeplusone.com/
# Or use R package ggseqlogo
```

## Plant-Specific Considerations

### 1. Large Gene Families

For families with 50+ members:

```bash
# Limit motifs to avoid redundancy
meme large_family.fasta \
  -mod anr \
  -nmotifs 15 \
  -minw 6 \
  -maxw 30 \
  -protein \
  -evt 0.01 \
  -cp 8
```

### 2. Domain-Specific Analysis

If interested in specific domains:

```bash
# Extract domain region first
# Example: Extract WRKY domain (~60aa)
grep -A 60 "WRKY" domain_annotation.txt > wrky_domains.fasta

meme wrky_domains.fasta \
  -mod anr \
  -nmotifs 5 \
  -minw 4 \
  -maxw 15 \
  -protein
```

### 3. Comparing Subgroups

For gene family subgroup analysis:

```bash
# Run MEME separately for each subgroup
# Group1
meme subgroup1.fasta -mod anr -nmotifs 10 -protein -oc group1_meme

# Group2
meme subgroup2.fasta -mod anr -nmotifs 10 -protein -oc group2_meme

# Compare motifs between groups
```

## Integration with Gene Family Pipeline

```bash
# Complete motif analysis pipeline

# 1. Extract protein sequences for family members
seqkit grep -f gene_ids.txt genome_proteins.fasta > family_proteins.fasta

# 2. Run MEME
meme family_proteins.fasta \
  -mod anr \
  -nmotifs 10 \
  -minw 6 \
  -maxw 50 \
  -protein \
  -oc meme_output \
  -cp 8

# 3. Annotate motifs with TOMTOM
tomtom -oc tomtom_output \
  -db $MEME_DB/motif_databases/Pfam.meme \
  meme_output/meme.txt

# 4. Visualize with TBtools or custom R script
```

## Output Interpretation

### E-value

Motif E-value indicates statistical significance:
- **E < 1e-10**: Highly significant, likely functional
- **E < 1e-3**: Significant, worth investigating
- **E > 0.1**: May be spurious

### Motif Width

Biological interpretation:
- **4-8 aa**: Short motif, possibly linear epitope
- **8-20 aa**: Typical protein domain core
- **20-50 aa**: Extended domain or structural element

### Sequence Logo

Height of letters = information content:
- Tall letters = conserved positions
- Short letters = variable positions

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| No motifs found | Sequences too diverse | Reduce `-nmotifs`, increase `-maxw` |
| Too many motifs | Parameters too relaxed | Increase `-evt` threshold |
| Motifs too short | Default `-minw` too low | Set `-minw 8` or higher |
| Motifs too long | `-maxw` too high | Set `-maxw 30` for proteins |
| Memory error | Too many sequences | Split into subgroups |

## Related Tools in MEME Suite

| Tool | Purpose |
|------|---------|
| **TOMTOM** | Compare motifs to databases |
| **FIMO** | Scan sequences for motif matches |
| **MAST** | Motif alignment and search |
| **DREME** | Short motif discovery |
| **MEME-ChIP** | Complete ChIP-seq motif analysis |

## References

- Bailey & Elkan (1994) Fitting a mixture model by expectation maximization to discover motifs in bipolymers
- MEME Suite: https://meme-suite.org/

## Related Tools

- **TBtools**: Plant-specific motif analysis and visualization
- **ggseqlogo**: R package for sequence logos
- **WebLogo**: Online sequence logo generator
