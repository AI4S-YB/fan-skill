# MCScanX (Multiple Collinearity Scan)

**Goal:** Detect syntenic blocks and collinear gene pairs within or between genomes
**Best for:** Genome duplication analysis, inter-species synteny comparison

## Prerequisites
- MCScanX installed (C++ tool) or JCVI MCScan (Python)
- BLASTP results (or DIAMOND for speed)
- GFF3 annotation file with gene positions

## Basic Usage

```bash
# 1. BLASTP all-vs-all between two species
blastp -query speciesA.pep -db speciesB.pep -out blast.txt -outfmt 6 -evalue 1e-5

# 2. MCScanX
MCScanX speciesA_speciesB_prefix
```

## Plant-Specific Notes
- For polyploid species (wheat, cotton, canola): run per-subgenome comparisons
- MCScanX outputs `.collinearity` (gene pairs) and `.html` (dot plot)
- Use `jcvi` Python version for multi-species synteny visualization
- Plant genomes have extensive synteny — expect >50% of genes in collinear blocks for close species

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "segmentation fault" | GFF format error | Check tab-delimited format |
| "no collinear blocks" | BLAST evalue too strict | Relax to 1e-3 |
