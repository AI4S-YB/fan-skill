# PAML: Phylogenetic Analysis by Maximum Likelihood

## Tool Overview

**Tool ID**: `paml`
**Category**: Molecular Evolution / Phylogenetics
**Purpose**: Maximum likelihood analysis of protein and DNA sequences, including Ka/Ks calculation
**Source**: http://abacus.gene.ucl.ac.uk/software/paml.html
**Version**: PAML 4.x

## When to Use

- Ka/Ks (dN/dS) ratio calculation for gene pairs
- Positive selection detection
- Molecular clock testing
- Ancestral sequence reconstruction

## Key Programs

| Program | Purpose |
|---------|---------|
| **yn00** | Pairwise Ka/Ks calculation (Yang-Nielsen 2000 method) |
| **codeml** | Codon-based analysis (site models, branch models) |
| **baseml** | Nucleotide substitution analysis |
| **basemlg** | Continuous gamma model for nucleotides |

## yn00: Pairwise Ka/Ks Calculation

### When to Use

- Gene family duplication analysis
- WGD gene pair evolution
- Small number of gene pairs (<100)

### Input Preparation

```bash
# yn00 requires codon-aligned sequences in PHYLIP format
# Step 1: Align CDS sequences by protein alignment
pal2nal.pl protein_alignment.fasta cds_sequences.fasta -output paml > aligned.phy

# Step 2: Create yn00 control file
cat > yn00.ctl << 'EOF'
seqfile = aligned.phy      * sequence file
outfile = yn00_results.txt * output file
verbose = 1                * output level
icode = 0                  * genetic code (0=universal)
weighting = 0              * weighting method
commonkappa = 0            * common kappa
commonf3x4 = 0             * common f3x4
END
EOF

# Step 3: Run yn00
yn00
```

### Output Interpretation

```
Sequence pair: Gene001 vs Gene002

Nei & Gojobori 1986. dN/dS = 0.2345
Yang & Nielsen 2000. dN/dS = 0.2567
  dN = 0.0234  dS = 0.0912

LWL85 method: dN/dS = 0.2456
```

**Key values**:
- **dN (Ka)**: Non-synonymous substitution rate
- **dS (Ks)**: Synonymous substitution rate
- **dN/dS (Ka/Ks)**: Selection indicator

### Ka/Ks Interpretation

| Ka/Ks | Interpretation | Biological Meaning |
|-------|----------------|-------------------|
| < 0.1 | Strong purifying selection | Highly conserved function |
| 0.1 - 0.5 | Moderate purifying selection | Functional constraint |
| 0.5 - 1.0 | Weak purifying selection | Relaxed constraint |
| ≈ 1.0 | Neutral evolution | No selection pressure |
| > 1.0 | Positive selection | Adaptive evolution |

## codeml: Advanced Analysis

### Site Models (Detecting Positive Selection Sites)

```bash
# Create control file for site model
cat > codeml.ctl << 'EOF'
seqfile = aligned.phy
treefile = tree.nwk
outfile = codeml_results.txt
runmode = 0        * user-defined tree
model = 0          * one ratio for all branches
NSsites = 3        * positive selection sites (M8 vs M7)
icode = 0          * universal genetic code
clock = 0          * no clock
fix_kappa = 0      * estimate kappa
kappa = 2          * initial kappa
fix_omega = 0      * estimate omega
omega = 0.5        * initial omega
getSE = 1          * get standard errors
END
EOF

codeml codeml.ctl
```

### Branch Models (Lineage-Specific Selection)

```bash
# Test for positive selection on specific branches
# model = 2 (two ratios), NSsites = 0
# Tree file must mark foreground branches with #1
# (Gene001#1, Gene002, Gene003);

codeml codeml_branch.ctl
```

## Plant-Specific Considerations

### 1. WGD Gene Pairs

```bash
# Extract WGD pairs from MCScanX output
awk '$5=="WGD" {print $2"\t"$3}' collinearity.out > wgd_pairs.txt

# For each pair, calculate Ka/Ks
while read gene1 gene2; do
  # Extract sequences
  seqkit grep -p "$gene1,$gene2" cds.fasta > pair.fasta
  
  # Align and run yn00
  # ... alignment steps ...
  yn00
done < wgd_pairs.txt
```

### 2. Multi-WGD Species

Plants with multiple WGD events show characteristic Ka/Ks peaks:

```r
# Plot Ka/Ks distribution to identify WGD peaks
kaks_data <- read.table("all_kaks_results.txt", header=TRUE)
hist(kaks_data$Ka_Ks, breaks=100, main="Ka/Ks Distribution",
     xlab="Ka/Ks", col="lightblue")

# Look for peaks:
# - Peak at Ka/Ks < 0.1: Recent WGD
# - Peak at Ka/Ks ~0.2-0.5: Medium-age WGD
# - Peak at Ka/Ks > 0.5: Ancient WGD
```

### 3. Genetic Code

```bash
# For plant chloroplast genes, use appropriate genetic code
# icode values in PAML:
# 0 = Universal (standard)
# 1 = Mammalian mitochondrial
# 2 = Yeast mitochondrial
# 3 = Mold mitochondrial
# 4 = Invertebrate mitochondrial
# 5 = Ciliate nuclear
# 6 = Echinoderm mitochondrial
# 9 = Echinoderm mitochondrial
# 10 = Alternative yeast nuclear
# 11 = Bacterial
# 12 = Alternative flatworm mitochondrial
# 13 = Chloroplast (for plant chloroplast genes!)
```

## Batch Processing

### For Large Gene Families (>100 pairs)

Use KaKs_Calculator instead of yn00 for efficiency:

```bash
# KaKs_Calculator is faster for large batches
# See kaks_calculator.md for details
```

### Automated Pipeline

```bash
#!/bin/bash
# batch_kaks.sh - Batch Ka/Ks calculation with PAML

PAIRS_FILE=$1
CDS_FILE=$2
PROT_FILE=$3

while read gene1 gene2; do
  # Extract sequences
  seqkit grep -p "$gene1" $PROT_FILE > prot1.fasta
  seqkit grep -p "$gene2" $PROT_FILE > prot2.fasta
  seqkit grep -p "$gene1" $CDS_FILE > cds1.fasta
  seqkit grep -p "$gene2" $CDS_FILE > cds2.fasta
  
  # Align proteins
  mafft --auto prot1.fasta prot2.fasta > aligned_prot.fasta
  
  # Pal2nal for codon alignment
  pal2nal.pl aligned_prot.fasta <(cat cds1.fasta cds2.fasta) \
    -output paml > aligned.phy
  
  # Run yn00
  echo "seqfile = aligned.phy" > yn00.ctl
  echo "outfile = ${gene1}_${gene2}_kaks.txt" >> yn00.ctl
  echo "verbose = 0, icode = 0, weighting = 0" >> yn00.ctl
  echo "commonkappa = 0, commonf3x4 = 0" >> yn00.ctl
  echo "END" >> yn00.ctl
  
  yn00
  
done < $PAIRS_FILE
```

## Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| "sequence error" | Sequences not codon-aligned | Use pal2nal for proper alignment |
| Ka/Ks = 999 | dS = 0 (no synonymous sites) | Pair too similar or alignment error |
| Ka/Ks = 0 | dN = 0 (identical proteins) | Check if sequences are identical |
| Negative Ka/Ks | Alignment issue | Re-align with proper codon awareness |
| No output | Control file syntax error | Check for spaces vs tabs, proper format |

## Integration with Gene Family Analysis

```bash
# Complete workflow: MCScanX → Ka/Ks analysis

# 1. Run MCScanX to identify WGD pairs
MCScanx gene_location

# 2. Extract WGD gene pairs
awk '/^Alignment/ {getline; for(i=1;i<=NF;i++) {if($i~/^[0-9]/) {split($i,a,"-"); print a[1], a[2]}}}' \
  gene_location.collinearity > wgd_pairs.txt

# 3. Calculate Ka/Ks for each pair
# Use batch script or KaKs_Calculator

# 4. Analyze results
# - Identify recent vs ancient WGD
# - Find genes under positive selection (Ka/Ks > 1)
# - Compare functional categories
```

## References

- Yang Z (2007) PAML 4: phylogenetic analysis by maximum likelihood
- Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution rates
- PAML Documentation: http://abacus.gene.ucl.ac.uk/software/paml.html

## Related Tools

- **KaKs_Calculator**: Faster batch Ka/Ks calculation
- **PAL2NAL**: Convert protein alignment to codon alignment
- **PAMLX**: GUI for PAML
- **HyPhy**: Alternative selection analysis
