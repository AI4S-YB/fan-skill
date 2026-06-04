# MCScanX: Synteny and Collinearity Analysis

## Tool Overview

**Tool ID**: `mcscanx`
**Category**: Comparative Genomics / Synteny Analysis
**Purpose**: Detect synteny blocks, whole-genome duplication (WGD) events, and gene duplication types
**Source**: https://github.com/wyp1125/MCScanX

## When to Use

- Gene family synteny analysis
- Whole-genome duplication detection
- Collinearity analysis between species
- Gene duplication type classification (WGD/tandem/proximal/dispersed)

## Key Parameters

| Parameter | Default | Recommended for Plants | Description |
|-----------|---------|------------------------|-------------|
| MATCH_SCORE | 50 | 50 | Score for a matching gene pair |
| GAP_PENALTY | -1 | **-25** | Penalty for introducing a gap |
| MIN_ALIGN | 5 | **5** | Minimum genes in a syntenic block |
| MAX_GAPS | 25 | 25 | Maximum gaps allowed in alignment |

### Critical Parameters for Plants

**GAP_PENALTY = -25** (more negative than default)

Why: Plant genomes have undergone multiple rounds of WGD and gene loss. A more negative gap penalty allows detection of syntenic blocks despite gene loss and rearrangement.

**MIN_ALIGN = 5**

Why: Minimum of 5 gene pairs provides statistical confidence while being sensitive enough to detect small syntenic blocks from ancient WGD.

## Input Files Required

### 1. Gene Location File (gene_location.txt)

Tab-separated format:
```
Chromosome	GeneID	Start	End
chr1	Gene001	1000	2000
chr1	Gene002	3000	4500
chr2	Gene003	1000	2000
```

Generate from GFF3:
```bash
awk '$3=="gene" {
  chr=$1;
  gsub(/ID=/, "", $9);
  gsub(/;.*/, "", $9);
  gene_id=$9;
  start=$4;
  end=$5;
  print chr"\t"gene_id"\t"start"\t"end
}' annotation.gff3 | sort -k1,1 -k3,3n > gene_location.txt
```

### 2. BLASTP Result File (gene_location.blast)

Self-BLAST for intra-genomic synteny, or cross-species BLAST for inter-genomic synteny:

```bash
# Make BLAST database
makeblastdb -in proteins.fasta -dbtype prot -out proteins_db

# Run BLASTP (all-vs-all)
blastp -query proteins.fasta \
  -db proteins_db \
  -outfmt 6 \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -out gene_location.blast
```

## Running MCScanX

```bash
# Basic syntax (input file prefix)
MCScanx gene_location

# This reads:
#   - gene_location.txt (gene locations)
#   - gene_location.blast (BLAST results)
# And produces:
#   - gene_location.collinearity (syntenic blocks)
#   - gene_location.tandem (tandem duplicates)
#   - gene_location.html (visualization)
```

## Output Interpretation

### 1. Collinearity File

```
############## Parameters ###############
MATCH_SCORE: 50
GAP_PENALTY: -25
MIN_ALIGN: 5
############## Statistics ###############
Total number of collinear blocks: 150
...
##########################################

Alignment 1: chr1 vs chr1
0-  152:  Gene001  Gene005  5.2e-80
0-  153:  Gene002  Gene006  2.1e-75
...
```

**Columns**:
- `0-`: Indicates whether genes are in same (0) or different orientation
- Gene pairs with E-value
- Syntenic block ID

### 2. Tandem File

Lists tandem duplicate genes:
```
Gene001	Gene002	2
Gene050	Gene051	Gene052	3
```
Number indicates tandem cluster size.

### 3. Duplication Type Classification

MCScanX classifies gene duplication types:

| Type | Description | Ka/Ks Pattern |
|------|-------------|---------------|
| **WGD** | Whole-genome duplication | Low-moderate Ka/Ks |
| **Tandem** | Adjacent genes | Variable Ka/Ks |
| **Proximal** | Nearby but not adjacent | Variable |
| **Dispersed** | Single copy or scattered | Variable |
| **Segmental** | Large duplicated blocks | Low-moderate Ka/Ks |

## Plant-Specific Workflow

### Detecting Ancient WGD

```bash
# Run MCScanX
MCScanx gene_location

# Extract WGD gene pairs
grep "WGD" gene_location.collinearity > wgd_pairs.txt

# Count WGD events
grep -c "WGD" gene_location.collinearity
```

### Cross-Species Synteny

For comparing two plant genomes:

```bash
# Combine gene locations from both species
cat species1_genes.txt species2_genes.txt > combined_genes.txt

# Cross-species BLASTP
blastp -query species1_proteins.fasta \
  -db species2_proteins.fasta \
  -outfmt 6 -evalue 1e-5 \
  -out cross_species.blast

# Run MCScanX
MCScanx combined_genes
```

### Visualization with Circos

```bash
# Generate Circos links from MCScanX output
awk 'NR>10 && /^Alignment/ {
  # Parse syntenic blocks for Circos
}' gene_location.collinearity > circos_links.txt
```

## Integration with Ka/Ks Analysis

```bash
# Extract syntenic gene pairs
awk '$5=="WGD" {print $2"\t"$3}' gene_location.collinearity > wgd_pairs.txt

# Calculate Ka/Ks for WGD pairs
# Use KaKs_Calculator or PAML
KaKs_Calculator -i wgd_pairs.axt -o kaks_results.txt
```

## Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| No syntenic blocks found | MIN_ALIGN too high | Reduce to 3-5 |
| Too many small blocks | GAP_PENALTY too lenient | Use -25 or lower |
| Missing known WGD | Gene annotation incomplete | Check BUSCO completeness |
| Memory error | Genome too large | Process by chromosome |

## Plant WGD Interpretation

### Expected Patterns

| WGD Type | Age | Ka/Ks Range | Pattern |
|----------|-----|-------------|---------|
| Recent WGD | <5 Mya | Ka/Ks < 0.1 | Many duplicate genes retained |
| Medium WGD | 5-50 Mya | Ka/Ks 0.1-0.5 | Partial gene loss |
| Ancient WGD | >50 Mya | Ka/Ks > 0.5 | Most duplicates lost, few retained |

### Multi-WGD Species

Plants often have multiple WGD events:
- **Brassica**: Triplication (~15 Mya) + paleohexaploidy
- **Maize**: Allotetraploidy (~12 Mya) + ancient WGD
- **Wheat**: Three rounds of polyploidy

Use `Ka/Ks distribution` to separate WGD peaks:
```r
# Plot Ka/Ks distribution to identify WGD peaks
hist(kaks_values, breaks=100, main="Ka/Ks Distribution")
# Look for peaks corresponding to different WGD events
```

## References

- Wang et al. (2012) MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity
- GitHub: https://github.com/wyp1125/MCScanX

## Related Tools

- **JCvi**: Python alternative for synteny analysis
- **SynVisio**: Interactive synteny visualization
- **Circos**: Whole-genome visualization including synteny
- **D-Genies**: Dot plot for genome comparisons
