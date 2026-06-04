# IQ-TREE: Phylogenetic Tree Construction

## Tool Overview

**Tool ID**: `iqtree`
**Category**: Phylogenetics
**Purpose**: Maximum-likelihood phylogenetic tree inference with automatic model selection
**Web**: http://www.iqtree.org/
**Version**: IQ-TREE 2 (recommended)

## When to Use

- Gene family phylogenetic analysis
- Species tree construction
- SNP-based population phylogeny
- Molecular evolution analysis

## Key Parameters

| Parameter | Default | Recommended | Description |
|-----------|---------|-------------|-------------|
| `-m` | MFP | **MFP** | Model selection (ModelFinder Plus) |
| `-B` | - | **1000** | Ultrafast bootstrap replicates |
| `--alrt` | - | **1000** | SH-aLRT test replicates |
| `-T` | 1 | **AUTO** | Number of threads |
| `--asc-corr` | - | **lewis** | ASC correction for SNP data |

### Critical Parameters for Plants

#### 1. Model Selection: `-m MFP`

IQ-TREE's ModelFinder Plus automatically selects the best substitution model:
- Tests 100+ models including mixture models
- For protein sequences: tests LG, WAG, JTT, etc. with rate heterogeneity
- For SNP data: **must add +ASC**

#### 2. ASC Correction for SNP Data

**Essential for plant SNP phylogenies!**

SNP data lacks invariant sites, causing biased likelihood estimates. Use:

```bash
# SNP data MUST use ASC correction
iqtree2 -s snp_alignment.fasta \
  -m MFP+ASC \
  --asc-corr lewis \
  -B 1000 --alrt 1000
```

**Why Lewis correction?**
- SNP data is a filtered version of full sequence data
- Invariant sites are excluded, violating model assumptions
- Lewis correction accounts for this ascertainment bias

#### 3. Bootstrap Support

Two support values recommended:

| Method | Threshold | Meaning |
|--------|-----------|---------|
| **UFBoot** (`-B`) | ≥ 95% | Strong support |
| **SH-aLRT** (`--alrt`) | ≥ 80% | Strong support |

Report both values: `Node support (UFBoot/SH-aLRT): 98/95`

## Command Examples

### Protein Sequences (Gene Family)

```bash
# Standard protein phylogeny
iqtree2 -s aligned_proteins.fasta \
  -m MFP \
  -B 1000 --alrt 1000 \
  -T 8 \
  --prefix gene_family

# Output files:
# .treefile - ML tree with support values
# .iqtree - Full report with model selection
# .log - Log file
# .ufboot - Bootstrap trees
```

### Coding Sequences (CDS)

```bash
# Codon-based phylogeny (more informative)
iqtree2 -s aligned_cds.fasta \
  -m MFP \
  -B 1000 --alrt 1000 \
  -st CODON \
  -T 8
```

### SNP Data (Population Phylogeny)

```bash
# SNP phylogeny with ASC correction
iqtree2 -s snp_alignment.fasta \
  -m MFP+ASC \
  --asc-corr lewis \
  -B 1000 --alrt 1000 \
  -T 8
```

### Partitioned Analysis (Multi-Gene)

```bash
# Create partition file (partitions.txt)
# [Gene1]
# gene1_region
# [Gene2]
# gene2_region

iqtree2 -s concatenated.fasta \
  -p partitions.txt \
  -m MFP+MERGE \
  -B 1000 --alrt 1000 \
  -T 16
```

## Model Selection Output

IQ-TREE reports the best-fit model:

```
Best-fit model according to BIC: LG+F+R5

Model components:
- LG: Le-Gascuel substitution matrix
-+F: Empirical amino acid frequencies
+R5: FreeRate model with 5 categories
```

**Common models for plants**:

| Data Type | Common Best-Fit Models |
|-----------|----------------------|
| Protein | LG+G4, LG+F+G4, WAG+G4 |
| Nucleotide | GTR+F+G4, TVM+F+G4 |
| Codon | GY+F+R3, MG+F+R3 |
| SNP | GTR+ASC+G4 |

## Plant-Specific Considerations

### 1. Polyploid Gene Families

In polyploids, distinguish:
- **Orthologs**: Same gene in different species (speciation)
- **Homeologs**: Duplicated genes from WGD (within species)
- **Paralogs**: General gene duplicates

**Strategy**: Include all copies, use tree topology to infer relationships

### 2. Large Gene Families

Plant gene families can have 50-200+ members:

```bash
# For large alignments, use faster approximations
iqtree2 -s large_family.fasta \
  -m MFP \
  -B 1000 \
  -nstop 500 \
  -T AUTO
```

### 3. Long Branches

Plant gene families may have highly divergent members:

```bash
# Use posterior mean site frequency (PMSF) model
iqtree2 -s alignment.fasta \
  -m LG+F+R5 \
  --pmsf model.params \
  -B 1000
```

## Output Interpretation

### Tree File Format

```
(GeneA:0.123,(GeneB:0.045,GeneC:0.067)90/85:0.089,GeneD:0.234);
```

Format: `(Node1:branch_length,Node2:branch_length)support_values`

### Support Value Interpretation

| UFBoot | SH-aLRT | Interpretation |
|--------|---------|----------------|
| ≥ 95 | ≥ 80 | Strong support - reliable clade |
| 70-94 | 60-79 | Moderate support - needs confirmation |
| < 70 | < 60 | Weak support - unreliable |

### Branch Length

Branch length = expected number of substitutions per site

- Long branch: Rapid evolution or divergence
- Short branch: Recent divergence or conserved sequence

## Visualization

### FigTree (Desktop)

```bash
# Open .treefile in FigTree GUI
figtree gene_family.treefile
```

### iTOL (Web)

1. Upload `.treefile` to https://itol.embl.de/
2. Customize colors, labels, and annotations
3. Export as publication-quality figure

### R (ggtree)

```r
library(ggtree)
tree <- read.tree("gene_family.treefile")
ggtree(tree) + geom_tiplab() + geom_nodelab()
```

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| Model selection fails | Too many invariant sites | Use `-mset` to limit models |
| Tree unresolved | Low phylogenetic signal | Add more loci or taxa |
| Long branch attraction | Divergent sequences | Add intermediate taxa |
| SNP tree wrong | Missing ASC correction | Add `-m MFP+ASC --asc-corr lewis` |
| Memory error | Large alignment | Use `-mem 32G` or reduce taxa |

## Integration with Gene Family Analysis

```bash
# Complete gene family workflow

# 1. Align sequences
mafft --auto --thread 8 proteins.fasta > aligned.fasta

# 2. Trim alignment
trimal -in aligned.fasta -out aligned_trimmed.fasta -automated1

# 3. Build tree
iqtree2 -s aligned_trimmed.fasta \
  -m MFP \
  -B 1000 --alrt 1000 \
  -T 8 \
  --prefix family_tree

# 4. Visualize and annotate
# Use iTOL with gene annotation file
```

## References

- Nguyen et al. (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies
- Minh et al. (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference
- Hoang et al. (2018) UFBoot2: Improving the ultrafast bootstrap approximation

## Related Tools

- **MAFFT**: Multiple sequence alignment
- **trimAl**: Alignment trimming
- **FigTree**: Tree visualization
- **iTOL**: Interactive tree visualization
- **RAxML-NG**: Alternative ML tree builder
