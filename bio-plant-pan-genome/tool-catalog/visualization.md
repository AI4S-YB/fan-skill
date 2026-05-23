# Pan-genome Visualization

**Goal:** Visualize pan-genome graphs, PAV distributions, core/variable gene composition
**Best for:** Publication-quality figures and exploratory analysis of plant pan-genomes

## ODGI Visualization

### Graph Layout and Rendering

```bash
# Build ODGI from GFA
odgi build -g pan_graph.gfa -o pan_graph.og

# Sort graph for better layout
odgi sort -i pan_graph.og -o pan_sorted.og -p s

# Generate 2D layout
odgi layout -i pan_sorted.og -o pan_layout.lay -t 16

# Generate PNG visualization
odgi draw -i pan_sorted.og \
  -c pan_layout.lay \
  -p pan_viz.png \
  -w 2047 -H 2047 \
  --show-strand

# Generate 1D path visualization (gene-level)
odgi extract -i pan_sorted.og \
  -r "gene:Os01g0100100" \
  -o gene_region.og

odgi viz -i gene_region.og \
  -o gene_region_viz.png \
  -x 2000 -y 200
```

## PAV Matrix Heatmap

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load PAV matrix
pav = pd.read_csv("pav_matrix.csv", index_col=0)

# Cluster by gene presence pattern
sns.clustermap(pav,
               method='ward',
               metric='jaccard',
               cmap='RdYlBu_r',
               figsize=(12, 16),
               xticklabels=True,
               yticklabels=False,
               cbar_kws={'label': 'Present (1) / Absent (0)'})
plt.title('Pan-genome PAV Landscape')
plt.savefig('pav_heatmap.png', dpi=300, bbox_inches='tight')
```

## Core/Variable Gene Composition

```python
import matplotlib.pyplot as plt

# Pie chart of gene categories
categories = ['Core\n(≥95%)', 'Soft-core\n(80-95%)',
              'Variable\n(15-80%)', 'Private\n(<15%)']
counts = [core_count, soft_core_count, variable_count, private_count]
colors = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c']

plt.figure(figsize=(8, 8))
plt.pie(counts, labels=categories, colors=colors,
        autopct='%1.1f%%', startangle=90,
        explode=(0, 0, 0, 0.1))
plt.title('Pan-genome Gene Composition')
plt.savefig('core_variable_pie.png', dpi=300, bbox_inches='tight')
```

## Genome-Specific PAV Profiles

```python
# Bar plot: number of private/variable genes per genome
import matplotlib.pyplot as plt
import pandas as pd

pav = pd.read_csv("pav_matrix.csv", index_col=0)
freq = pav.sum(axis=1) / pav.shape[1]

def classify(f):
    if f >= 0.95: return 'Core'
    elif f >= 0.80: return 'Soft-core'
    elif f >= 0.15: return 'Variable'
    else: return 'Private'

gene_classes = freq.apply(classify)

# Stacked bar per genome
profile = pd.DataFrame(index=pav.columns)
for cat in ['Core', 'Soft-core', 'Variable', 'Private']:
    cat_genes = gene_classes[gene_classes == cat].index
    profile[cat] = pav.loc[cat_genes].sum()

profile.plot(kind='bar', stacked=True, figsize=(12, 6),
             color=['#2ecc71', '#3498db', '#f39c12', '#e74c3c'])
plt.title('Gene Composition per Genome')
plt.ylabel('Number of Genes')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('genome_composition.png', dpi=300)
```

## Plant-Specific Figure Notes

- Multi-panel figure recommended: (a) ODGI graph visualization, (b) PAV heatmap, (c) core/variable pie chart
- For polyploids: color-code nodes by subgenome origin
- For publication: use Circos plots to show PAV distribution across chromosomes
- Include BUSCO completeness bar per genome to distinguish assembly issues from true PAV

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| ODGI draw crashes | Graph too large | Extract chromosome-level subgraphs |
| Heatmap too dense | Too many genes | Sample 5000 random genes or filter to variable only |
| Pie chart misleading | Private genes inflated by one poor genome | Remove outlier genome first |
| Circos labels overlap | Too many chromosomes | Group into chromosome-level tracks |
