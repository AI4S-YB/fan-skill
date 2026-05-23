# GRN Visualization — Network Plots for Plant Regulatory Networks

**Goal:** Publication-quality visualizations of gene regulatory networks:
TF-target subnetworks, hub gene ranking, module structure, and condition comparisons.

**Best for:** All GRN analysis results.

## Prerequisites

- R 4.0+: igraph, ggraph, ggplot2, tidygraph, RColorBrewer
- Python 3.8+: networkx, matplotlib, pyvis (interactive HTML), plotly
- Input: edge list, hub gene list, module assignments, centrality scores

## TF-Target Subnetwork Plot (Python + pyvis)

```python
import pandas as pd
import networkx as nx
from pyvis.network import Network

# Load edges and hub genes
edges = pd.read_csv("grn_edges_genie3.tsv", sep="\t")
hubs = pd.read_csv("hub_genes.tsv", sep="\t")
hub_genes = set(hubs[hubs["is_hub"]]["gene"])

# Filter to top 500 edges for visualization
top_edges = edges.nlargest(500, "importance")

# Build subgraph
G = nx.from_pandas_edgelist(
    top_edges,
    source="TF",
    target="target",
    edge_attr="importance",
    create_using=nx.DiGraph()
)

# Interactive network
net = Network(height="800px", width="100%",
              directed=True, notebook=False)

for node in G.nodes():
    if node in hub_genes:
        # Hub genes: large red nodes
        net.add_node(node, label=node, size=25, color="#e74c3c",
                     title=f"Hub: {node}")
    elif node in tf_list:
        # Other TFs: medium blue nodes
        net.add_node(node, label=node, size=12, color="#3498db",
                     title=f"TF: {node}")
    else:
        # Target genes: small grey nodes
        net.add_node(node, label=node, size=5, color="#95a5a6",
                     title=f"Target: {node}")

for u, v, d in G.edges(data=True):
    net.add_edge(u, v, value=d.get("importance", 1))

net.show("grn_subnetwork.html")
```

## Hub Gene Ranking Plot (Python)

```python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

hubs = pd.read_csv("hub_genes.tsv", sep="\t")
tf_hubs = hubs[hubs["is_tf"]].nlargest(20, "hub_score")

fig, ax = plt.subplots(figsize=(10, 6))

colors = plt.cm.viridis(np.linspace(0, 1, len(tf_hubs)))
ax.barh(range(len(tf_hubs)), tf_hubs["hub_score"].values,
        color=colors[::-1])
ax.set_yticks(range(len(tf_hubs)))
ax.set_yticklabels(tf_hubs["gene"].values)
ax.set_xlabel("Hub Score (combined Z-score)")
ax.set_title("Top 20 Hub Transcription Factors")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig("hub_tf_ranking.pdf", dpi=150, bbox_inches="tight")
```

## Module Structure Visualization (R)

```r
library(igraph)
library(ggraph)
library(ggplot2)

# Load edges and module assignments
edges <- read.table("grn_edges_genie3.tsv", header = TRUE, sep = "\t")
modules <- read.table("gene_modules.tsv", header = TRUE, sep = "\t")

# Build graph for a specific module
mod_genes <- modules$gene[modules$module == 0]
mod_edges <- edges[edges$TF %in% mod_genes & edges$target %in% mod_genes, ]
g <- graph_from_data_frame(mod_edges, directed = TRUE)

ggraph(g, layout = "fr") +
  geom_edge_link(aes(alpha = importance), color = "grey60") +
  geom_node_point(aes(color = name %in% tf_list), size = 3) +
  scale_color_manual(values = c("grey70", "#e74c3c"),
                     labels = c("Target", "TF")) +
  labs(title = sprintf("Module 0 (%d genes, %d edges)",
       length(mod_genes), nrow(mod_edges))) +
  theme_void() +
  theme(legend.position = "bottom")

ggsave("module_0_network.pdf", width = 10, height = 10)
```

## Condition Comparison Heatmap (Python)

```python
# Compare regulon activity across conditions
aucell = pd.read_csv("aucell_output.tsv", sep="\t", index_col=0)

# Z-score by row (regulon)
aucell_z = aucell.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Top regulons by variance
top_regs = aucell_z.var(axis=1).nlargest(20).index

fig, ax = plt.subplots(figsize=(12, 8))
im = ax.imshow(aucell_z.loc[top_regs], aspect="auto", cmap="RdBu_r",
               vmin=-2, vmax=2)
plt.colorbar(im, label="Regulon Activity (Z-score)")
ax.set_xticks(range(aucell_z.shape[1]))
ax.set_xticklabels(aucell_z.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(top_regs)))
ax.set_yticklabels(top_regs, fontsize=8)
ax.set_title("Top 20 Regulons Activity by Condition")
plt.tight_layout()
plt.savefig("regulon_activity_heatmap.pdf", dpi=150, bbox_inches="tight")
```

## Plant-Specific Visualization Notes

- **Chromosome coordinates**: If TF and target genomic positions are known,
  use Circos plots to overlay GRN edges on chromosome ideograms —
  especially informative for polyploids with distinct subgenomes.

- **TF family color coding**: Assign consistent colors to major plant TF families:
  MYB=red, WRKY=blue, NAC=green, bHLH=orange, AP2/ERF=purple, others=grey.

- **Condition comparison**: Plant GRN visualizations most impactful when showing
  network rewiring between conditions (e.g., control vs. drought, or tissue A vs. B).

- **Interactive exploration**: Use pyvis HTML exports for exploratory analysis.
  Include TF family, centrality, and expression level as node attributes for
  hover tooltips.
