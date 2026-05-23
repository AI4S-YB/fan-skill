# PPI Network Analysis (STRING / Cytoscape)

**Goal:** Map differentially abundant proteins to protein-protein interaction networks
**Best for:** Understanding functional relationships among DEPs in plant systems

## Prerequisites

- STRING database (https://string-db.org/)
- Cytoscape (https://cytoscape.org/) — optional for advanced visualization
- List of differentially abundant proteins with significance values

## STRING DB Analysis

### Web Interface

1. Go to https://string-db.org/
2. Select organism: *Arabidopsis thaliana*, *Oryza sativa*, *Zea mays*, etc.
3. Input protein list (UniProt IDs or locus tags)
4. Configure:
   - Network type: "full STRING network"
   - Required score: 0.4 (medium confidence) or 0.7 (high confidence)
   - Max interactors: 0 (no limit) for 1st shell, or 10-50 for 2nd shell
5. Export as TSV or Cytoscape-compatible format

### API-Based (R)

```r
library(STRINGdb)

# Initialize STRING database for Arabidopsis
string_db <- STRINGdb$new(
  version = "11.5",
  species = 3702,  # Arabidopsis thaliana
  score_threshold = 400,
  input_directory = ""
)

# Map protein IDs
dep_mapped <- string_db$map(deg_results, "protein_id",
                            removeUnmappedRows = TRUE)

# Get network
string_db$plot_network(dep_mapped$STRING_id)
```

## Cytoscape Network Analysis

### Import STRING Network

1. File > Import > Network from Public Databases
2. Select STRING: protein query
3. Enter protein IDs
4. Select species and confidence cutoff
5. Download and import

### Network Clustering (MCODE)

```
Apps > MCODE > Analyze Network

Parameters:
- Degree cutoff: 2
- Node score cutoff: 0.2
- K-core: 2
- Max depth: 100

Typical result: 3-10 clusters representing functional modules
```

### Functional Enrichment on Network

```
Apps > ClueGO or BiNGO

Parameters for ClueGO:
- Ontology: GO Biological Process (Plant)
- Evidence: all
- Statistical test: right-sided hypergeometric
- Correction: Bonferroni step down
- Min genes: 3, Min percentage: 4%
```

## Plant-Specific Network Considerations

### Available STRING Plant Species

| Species | STRING Taxon ID | Protein Count |
|---------|----------------|---------------|
| Arabidopsis thaliana | 3702 | 27,655 |
| Oryza sativa japonica | 39947 | 35,000+ |
| Zea mays | 4577 | 39,000+ |
| Glycine max | 3847 | 46,000+ |

### Non-Model Plant Strategy

```python
# Map proteins to STRING via orthology
import pandas as pd
from Bio import SeqIO

# 1. BLAST non-model proteins to model species
# blastp -query non_model.faa -db arabidopsis -out blast.tsv

# 2. Map based on reciprocal best hits
blast_results = pd.read_csv("blast.tsv", sep="\t", header=None)
blast_results.columns = ["query", "target", "identity", "length",
                          "mismatch", "gap", "qstart", "qend",
                          "tstart", "tend", "evalue", "score"]

# Filter: identity > 40%, coverage > 50%
best_hits = blast_results[
    (blast_results["identity"] > 40) &
    (blast_results["length"] / blast_results.groupby("query")["length"].transform("max") > 0.5)
]
```

## Network Metrics to Report

```python
import networkx as nx

# Load network
G = nx.read_edgelist("string_network.txt")

# Metrics
print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")
print(f"Clustering coefficient: {nx.average_clustering(G):.3f}")
print(f"Network density: {nx.density(G):.4f}")

# Hub proteins (top 10 by degree)
degrees = dict(G.degree())
hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:10]
print("Top hubs:", hubs)

# Betweenness centrality
bc = nx.betweenness_centrality(G)
top_bc = sorted(bc.items(), key=lambda x: x[1], reverse=True)[:10]
print("Top betweenness:", top_bc)
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No interactions found" | Protein IDs not recognized by STRING | Convert to STRING-recognized IDs |
| Network too sparse | Too few input proteins or high confidence | Lower confidence score to 0.4 |
| Network too dense (>1000 edges) | Too many proteins or low confidence | Filter to significant proteins only, raise score |
| Non-model species: poor mapping | Too-divergent proteins | Use OrthoFinder-based mapping |
| Cytoscape freezes | Too many nodes/edges | Filter network to top 200 nodes by degree |
