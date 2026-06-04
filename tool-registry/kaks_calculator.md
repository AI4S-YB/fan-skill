# KaKs_Calculator: Batch Ka/Ks Calculation

## Tool Overview

**Tool ID**: `kaks_calculator`
**Category**: Molecular Evolution
**Purpose**: Calculate Ka (non-synonymous) and Ks (synonymous) substitution rates for gene pairs
**Source**: https://sourceforge.net/projects/kakscalculator2/
**Version**: KaKs_Calculator 2.0 (recommended)

## When to Use

- Large-scale Ka/Ks calculation (>100 gene pairs)
- WGD gene pair analysis
- Gene family evolution studies
- Faster alternative to PAML for batch processing

## Key Parameters

| Parameter | Options | Recommended | Description |
|-----------|---------|-------------|-------------|
| `-m` | MLWL, LPB, LWL, NG, GY, etc. | **MLWL** | Calculation method |
| `-i` | - | - | Input AXT format file |
| `-o` | - | - | Output file |
| `-c` | 0-16 | 0 (Universal) | Genetic code |

### Calculation Methods

| Method | Full Name | Description | Recommended For |
|--------|-----------|-------------|-----------------|
| **MLWL** | Modified Li-Wu-Luo | Maximum likelihood version | **Plant gene families** |
| **LPB** | Li et al. | Probability-based | General use |
| **LWL** | Li-Wu-Luo | Classic method | Quick estimation |
| **NG** | Nei-Gojobori | Simple counting | Approximation |
| **GY** | Goldman-Yang | Codon-based ML | More accurate, slower |

## Input Format: AXT

KaKs_Calculator requires AXT format input:

```
>Gene001-Gene002
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
>Gene003-Gene004
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
```

**Format specification**:
- Header line: `>GenePairName` (no spaces)
- First sequence line: Gene 1 CDS
- Second sequence line: Gene 2 CDS
- Sequences must be codon-aligned (same length, multiple of 3)

## Command Examples

### Basic Usage

```bash
# Calculate Ka/Ks for all gene pairs
KaKs_Calculator -i gene_pairs.axt -o kaks_results.txt -m MLWL
```

### Batch Processing from MCScanX Output

```bash
#!/bin/bash
# Step 1: Extract gene pairs from MCScanX collinearity file
awk '/^[0-9]/ {
  # Parse gene pairs from collinearity output
  gene1 = $2;
  gene2 = $3;
  print gene1"\t"gene2
}' gene_location.collinearity | sort -u > gene_pairs.txt

# Step 2: Extract CDS sequences for each pair
while read gene1 gene2; do
  echo ">$gene1-$gene2"
  seqkit grep -p "$gene1" cds.fasta | tail -n +2 | tr -d '\n'
  echo ""
  seqkit grep -p "$gene2" cds.fasta | tail -n +2 | tr -d '\n'
  echo ""
done < gene_pairs.txt > gene_pairs.axt

# Step 3: Run KaKs_Calculator
KaKs_Calculator -i gene_pairs.axt -o kaks_results.txt -m MLWL
```

### With Codon Alignment (Recommended)

For accurate Ka/Ks, use codon-aware alignment:

```bash
#!/bin/bash
# Complete pipeline with codon alignment

PAIRS_FILE=$1
CDS_FILE=$2
PROT_FILE=$3
OUTPUT_DIR=$4

mkdir -p $OUTPUT_DIR
> all_pairs.axt

while read gene1 gene2; do
  # Extract protein sequences
  seqkit grep -p "$gene1" $PROT_FILE > temp_prot1.fasta
  seqkit grep -p "$gene2" $PROT_FILE > temp_prot2.fasta
  
  # Extract CDS sequences
  seqkit grep -p "$gene1" $CDS_FILE > temp_cds1.fasta
  seqkit grep -p "$gene2" $CDS_FILE > temp_cds2.fasta
  
  # Align proteins
  cat temp_prot1.fasta temp_prot2.fasta | \
    mafft --auto - > aligned_prot.fasta 2>/dev/null
  
  # Convert to codon alignment using PAL2NAL
  cat temp_cds1.fasta temp_cds2.fasta > temp_cds.fasta
  pal2nal.pl aligned_prot.fasta temp_cds.fasta -output axt >> all_pairs.axt
  
done < $PAIRS_FILE

# Run KaKs_Calculator
KaKs_Calculator -i all_pairs.axt -o $OUTPUT_DIR/kaks_results.txt -m MLWL
```

## Output Format

```
Sequence     Ka      Ks     Ka/Ks   P-Value(Fisher)  Length
Gene001-002  0.0234  0.0912 0.2567  0.0012           1500
Gene003-004  0.1567  0.0834 1.8780  0.0234           1200
```

**Key columns**:
- **Ka**: Non-synonymous substitution rate
- **Ks**: Synonymous substitution rate
- **Ka/Ks**: Selection indicator (dN/dS ratio)
- **P-Value(Fisher)**: Statistical significance

## Ka/Ks Interpretation

| Ka/Ks | Selection Type | Biological Meaning |
|-------|----------------|-------------------|
| < 0.1 | Strong purifying | Highly conserved function |
| 0.1 - 0.3 | Purifying selection | Functional constraint |
| 0.3 - 0.5 | Weak purifying | Some constraint |
| 0.5 - 1.0 | Relaxed selection | Functional drift |
| ≈ 1.0 | Neutral evolution | No selection |
| > 1.0 | Positive selection | Adaptive evolution |
| > 2.0 | Strong positive selection | Rapid adaptation |

## Plant-Specific Analysis

### WGD Age Estimation

```r
# Ks values correlate with WGD timing
# Approximate formula: T = Ks / (2r)
# where r = substitution rate per site per year

# For plants, r ≈ 6.5e-9 (nuclear genes)
# T = Ks / (2 × 6.5e-9) = Ks / 1.3e-8

# Example: Ks peak at 0.2
T = 0.2 / 1.3e-8  # ≈ 15.4 million years

# Plot Ks distribution to identify WGD peaks
kaks <- read.table("kaks_results.txt", header=TRUE)
hist(kaks$Ks, breaks=100, main="Ks Distribution",
     xlab="Ks", col="lightblue", xlim=c(0,2))
```

### Multi-WGD Species Analysis

Plants often have multiple WGD events with distinct Ks peaks:

```r
# Density plot to identify multiple WGD peaks
kaks <- read.table("kaks_results.txt", header=TRUE)
plot(density(kaks$Ks, na.rm=TRUE), main="Ks Density Plot",
     xlab="Ks", col="blue", lwd=2)

# Look for peaks at different Ks values:
# Ks ~ 0.05-0.1: Recent WGD (<10 Mya)
# Ks ~ 0.2-0.5: Medium-age WGD (15-40 Mya)  
# Ks ~ 0.6-1.5: Ancient WGD (>50 Mya)
```

### Subgenome-Specific Analysis

For polyploids (e.g., wheat, cotton):

```bash
# Separate homeolog pairs by subgenome
grep "A_genome-D_genome" gene_pairs.txt > AD_pairs.txt
grep "A_genome-B_genome" gene_pairs.txt > AB_pairs.txt

# Calculate Ka/Ks separately
KaKs_Calculator -i AD_pairs.axt -o AD_kaks.txt -m MLWL
KaKs_Calculator -i AB_pairs.axt -o AB_kaks.txt -m MLWL

# Compare Ka/Ks distributions
# Different distributions may indicate subgenome dominance
```

## Integration with Gene Family Workflow

```bash
# Complete workflow: Gene Family → Synteny → Ka/Ks

# 1. Identify gene family members
hmmsearch --tblout family.tbl domain.hmm proteins.fasta
cut -f1 family.tbl | grep -v "^#" | sort -u > family_ids.txt

# 2. Extract family proteins and run MCScanX
seqkit grep -f family_ids.txt proteins.fasta > family_proteins.fasta
# ... MCScanX workflow ...

# 3. Get syntenic pairs within family
grep "WGD\|segmental" collinearity.out | awk '{print $2"\t"$3}' > family_pairs.txt

# 4. Calculate Ka/Ks
# ... batch Ka/Ks calculation ...

# 5. Identify rapidly evolving genes
awk '$4 > 1.0' kaks_results.txt > positive_selection_genes.txt
```

## Performance Comparison

| Tool | Gene Pairs | Speed | Accuracy | Best For |
|------|-----------|-------|----------|----------|
| **KaKs_Calculator** | >100 | Fast | Good | Large batches |
| **PAML yn00** | <100 | Medium | High | Small sets, publication |
| **PAML codeml** | Variable | Slow | Highest | Detailed analysis |

## Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| Ka/Ks = 999 | Ks = 0 (identical or very similar) | Check sequence identity |
| Ka/Ks = 0 | Ka = 0 (no non-syn changes) | Normal for conserved genes |
| Negative values | Alignment or frame error | Re-align with codon awareness |
| "AXT format error" | Wrong input format | Check header and sequence lines |
| Empty output | No valid pairs | Check CDS extraction |

## References

- Wang et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating gamma series methods
- Zhang et al. (2006) KaKs_Calculator: calculating Ka and Ks through model selection

## Related Tools

- **PAML**: Alternative with codeml for advanced models
- **PAL2NAL**: Convert protein alignment to codon alignment
- **KaKs_plot**: Visualization of Ka/Ks distributions
