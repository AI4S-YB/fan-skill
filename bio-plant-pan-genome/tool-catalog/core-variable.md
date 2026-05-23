# Core/Variable Gene Classification

**Goal:** Classify genes into core, soft-core, variable, and private categories based on presence/absence across the pan-genome
**Best for:** Understanding genome composition and identifying dispensable/essential gene sets

## Prerequisites

- Gene presence/absence matrix (genes x genomes, binary or count)
- Orthogroup or gene family clustering results
- Python 3.6+ with pandas, numpy

## Basic Usage

```python
import pandas as pd
import numpy as np

# Load presence/absence matrix
pav = pd.read_csv("pav_matrix.csv", index_col=0)  # genes x genomes

# Calculate presence frequency
freq = pav.sum(axis=1) / pav.shape[1]

# Classify
def classify_core(frequency):
    if frequency >= 0.95:
        return "core"
    elif frequency >= 0.80:
        return "soft_core"
    elif frequency >= 0.15:
        return "variable"
    else:
        return "private"

gene_class = freq.apply(classify_core)

# Summary
print(gene_class.value_counts())
```

## Classification From Orthogroups

```bash
# Using OrthoFinder output
python3 << 'EOF'
import pandas as pd

# OrthoFinder Orthogroups.GeneCount.tsv
og = pd.read_csv("Orthogroups.GeneCount.tsv", sep="\t", index_col=0)
presence = (og > 0).astype(int)
freq = presence.sum(axis=1) / presence.shape[1]

for cat, thresh in [("core", 0.95), ("soft_core", 0.80),
                    ("variable", 0.15)]:
    count = (freq >= thresh).sum()
    print(f"{cat}: {count} orthogroups")
EOF
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| core threshold | >= 95% of genomes | Standard pan-genome definition |
| soft-core threshold | 80-95% | Genes present in most but not all |
| variable threshold | 15-80% | Dispensable genome |
| private threshold | < 15% | Rare or strain-specific genes |

## Plant-Specific Notes

- Plant NBS-LRR resistance genes are predominantly "variable" — expect 30-50% PAV
- Secondary metabolite gene clusters often classified as "variable" or "private"
- If >15% of genes are "private", check for annotation incompleteness in some genomes
- Polyploids: classify per subgenome, then compare patterns between subgenomes

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Core genome too small (<20%) | Poor orthogroup clustering, or too few genomes | Relax clustering parameters, add more genomes |
| Core genome too large (>80%) | Too few diverse genomes included | Add wild relatives or divergent accessions |
| Private genes = 0% | All genomes are too similar, or annotation pipeline identical | Include more diverse accessions |
| Inflation of private genes | One incomplete genome assembly | Check BUSCO completeness per genome |
