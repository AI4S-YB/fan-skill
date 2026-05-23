# Network Analysis — Topology, Centrality, and Community Detection

**Goal:** Analyze inferred GRN topology — identify hub nodes, detect functional
modules, and quantify network properties.

**Best for:** Any inferred GRN (GENIE3, SCENIC, or WGCNA output).

## Prerequisites

- Python 3.8+: networkx, pandas, numpy, scipy, python-igraph (or leidenalg)
- R 4.0+: igraph, WGCNA, dynamicTreeCut
- Input: edge list (source, target, weight) from any network inference method

## Centrality Analysis — Hub Gene Detection

```python
import pandas as pd
import networkx as nx
import numpy as np

# Load edge list
edges = pd.read_csv("grn_edges_genie3.tsv", sep="\t")
# Columns: TF (source), target, importance (weight)

# Build directed network
G = nx.from_pandas_edgelist(
    edges,
    source="TF",
    target="target",
    edge_attr="importance",
    create_using=nx.DiGraph()
)

print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")

# Compute centrality metrics
degree_cent = nx.degree_centrality(G)
between_cent = nx.betweenness_centrality(G, k=min(1000, G.number_of_nodes()))
eigen_cent = nx.eigenvector_centrality_numpy(G, weight="importance")

# Combine metrics into DataFrame
centrality_df = pd.DataFrame({
    "gene": list(degree_cent.keys()),
    "degree_centrality": list(degree_cent.values()),
    "betweenness_centrality": list(between_cent.values()),
    "eigenvector_centrality": list(eigen_cent.values()),
})

# Identify hub genes: top 5% in both degree and betweenness
p95_deg = centrality_df["degree_centrality"].quantile(0.95)
p95_bet = centrality_df["betweenness_centrality"].quantile(0.95)

hubs = centrality_df[
    (centrality_df["degree_centrality"] >= p95_deg) &
    (centrality_df["betweenness_centrality"] >= p95_bet)
].sort_values("degree_centrality", ascending=False)

print(f"Hub genes: {len(hubs)}")
hubs.to_csv("hub_genes.tsv", sep="\t", index=False)

# TF-specific hubs
tf_hubs = hubs[hubs["gene"].isin(tf_list)]
print(f"Hub TFs: {len(tf_hubs)}")
```

## Community Detection — Functional Modules

```python
import community as community_louvain  # python-louvain
# Or: from networkx.algorithms.community import greedy_modularity_communities

# Convert to undirected for community detection
G_undirected = G.to_undirected()

# Louvain community detection
partition = community_louvain.best_partition(G_undirected, weight="importance")

# Assign module membership
module_df = pd.DataFrame(
    partition.items(), columns=["gene", "module"]
)
module_df["module_size"] = module_df.groupby("module")["gene"].transform("count")

print(f"Modules: {module_df['module'].nunique()}")
print(f"Module sizes:\n{module_df['module_size'].describe()}")

# Filter to meaningful modules (>= 10 genes)
valid_modules = module_df[module_df["module_size"] >= 10]
module_df.to_csv("gene_modules.tsv", sep="\t", index=False)

# Module-TF enrichment: which TFs regulate each module?
for mod, mod_genes in module_df.groupby("module")["gene"]:
    mod_edges = edges[edges["target"].isin(mod_genes)]
    top_tfs = mod_edges.groupby("TF")["importance"].sum().nlargest(5)
    print(f"Module {mod} ({len(mod_genes)} genes) — top TFs: "
          f"{', '.join(top_tfs.index)}")
```

## Network-Level Statistics

```python
# Global network properties
print(f"Density: {nx.density(G):.6f}")
print(f"Average clustering: {nx.average_clustering(G.to_undirected()):.4f}")
print(f"Average shortest path: {nx.average_shortest_path_length(G.to_undirected()):.2f}")

# Out-degree distribution (TFs regulate how many targets?)
out_degrees = [d for n, d in G.out_degree()]
print(f"Mean out-degree: {np.mean(out_degrees):.1f}")
print(f"Median out-degree: {np.median(out_degrees):.1f}")
print(f"Max out-degree: {np.max(out_degrees)}")
```

## Plant-Specific Notes

- **Hub TFs in plants**: Expect top hub TFs to belong to plant-specific families
  (MYB, WRKY, NAC). If no plant-specific TFs appear in the top 20 hubs, check TF
  annotation completeness or network quality.

- **Module interpretation**: Plant co-expression modules often correspond to:
  photosynthesis, secondary metabolism, hormone signaling, stress response,
  and developmental processes. Cross-reference with GO enrichment.

- **Polyploid module structure**: In allopolyploids, functional modules may contain
  homeologous gene pairs from different subgenomes. If A-genome and B-genome
  paralogs consistently fall in different modules, this suggests subfunctionalization.

- **Drought/stress networks**: Under stress conditions, network rewiring is common.
  Hub genes may shift (e.g., from MYB in control to WRKY in drought). Compare
  hub lists and module assignments across conditions.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Graph has no edges" | Edge threshold too strict | Lower importance threshold |
| Memory error on betweenness | Network too large | Use k-approximation (k=1000) |
| All genes in one module | Resolution parameter too low | Increase Louvain resolution; try Leiden algorithm |
| Hub genes are all one family | TF list biased or artifact | Check if one TF family dominates the TF list itself |
| Negative eigenvector centrality | Disconnected graph | Check for isolated nodes; use largest component |
