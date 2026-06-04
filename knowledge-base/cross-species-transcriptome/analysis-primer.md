# Cross-Species Transcriptome Analysis Primer

## Analysis Overview

Cross-species transcriptome comparison enables understanding of gene expression evolution, tissue-specificity conservation, and regulatory divergence across species. This analysis framework integrates ortholog mapping, expression normalization, batch correction, and comparative statistical methods.

## Key Analysis Scenarios

### 1. Expression Conservation Analysis
- Identify genes with conserved expression patterns across species
- Quantify expression divergence using correlation metrics
- Detect tissue-specific expression shifts

### 2. Regulatory Evolution
- Link sequence evolution (dN/dS) with expression divergence
- Identify genes under expression constraint
- Detect positively selected expression patterns

### 3. Tissue-Specificity Comparison
- Calculate tissue-specificity indices (tau, SPM)
- Compare specificity conservation across orthologs
- Identify tissue-specific gene gains/losses

## Input Requirements

### Essential Data
- **Transcriptome data**: Expression matrices (TPM/FPKM/counts) for multiple species
- **Ortholog mapping**: One-to-one ortholog relationships between species
- **Species phylogeny**: Divergence times and evolutionary relationships

### Recommended Data
- **Tissue metadata**: Sample-to-tissue annotations
- **Codon alignments**: For dN/dS calculation (optional)
- **Gene annotations**: Functional categories for enrichment

## Critical Considerations

### Batch Effects (CRITICAL)
Different species data often come from different labs/platforms. **Always apply batch correction** using ComBat or similar methods, with species as the batch variable.

### Expression Normalization
- Use TPM for cross-species comparisons (accounts for gene length differences)
- Apply quantile normalization if distribution differences are severe
- Log-transform after adding pseudocount

### Ortholog Mapping Quality
- Prefer OrthoFinder for de novo ortholog identification
- Use biomaRt for well-annotated species in Ensembl
- Apply reciprocal best hit (RBH) as a simple fallback

## Expected Outputs

| Output | Description |
|--------|-------------|
| Ortholog expression matrix | Normalized expression for 1:1 orthologs |
| Tissue-specificity scores | tau/SPM indices per gene per species |
| Expression conservation | Cross-species correlation coefficients |
| Divergence estimates | Expression evolutionary rates |
| Visualization | PCA/t-SNE plots, heatmaps, scatter plots |

## Tools Overview

| Category | Primary Tools | Alternatives |
|----------|---------------|--------------|
| Ortholog identification | OrthoFinder | biomaRt, OrthoDB |
| Expression quantification | STAR + Salmon | HISAT2 + featureCounts |
| Batch correction | ComBat (sva) | MNN (batchelor) |
| Dimensionality reduction | PCA, t-SNE, UMAP | PCoA |
| Tissue specificity | tspex (tau, SPM) | Custom R scripts |
| Differential expression | limma | DESeq2, edgeR |
| Phylogenetic analysis | phytools, ape | PAML |

## References

1. Brawand D, et al. (2011) The evolution of gene expression levels in mammalian organs. Nature 478:343-348.
2. Chen J, et al. (2012) Tissue-specific comparison of gene expression in human and mouse. PLoS One 7:e34053.
3. Kryuchkova-Mostacci N, Robinson-Rechavi M (2017) A benchmark of gene expression tissue-specificity metrics. Brief Bioinform 18:205-214.
