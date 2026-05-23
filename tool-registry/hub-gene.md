# Hub Gene Analysis — Identify and Characterize Regulatory Hubs

**Goal:** Identify master regulator TFs (hub genes) from the GRN by integrating
centrality metrics, expression patterns, and functional annotations.

**Best for:** After GENIE3/SCENIC network inference — prioritize TFs for
experimental validation.

## Prerequisites

- Python 3.8+: pandas, numpy, networkx, scipy.stats
- Input: edge list + centrality results from `network-analysis.md`
- TF family annotation from `tf-database.md`
- Expression matrix (for hub-target correlation validation)

## Hub Gene Identification Pipeline

```python
import pandas as pd
import numpy as np
from scipy.stats import zscore

# Load centrality results
cent = pd.read_csv("centrality_scores.tsv", sep="\t")

# Z-score normalization for each metric
cent["deg_zscore"] = zscore(cent["degree_centrality"])
cent["bet_zscore"] = zscore(cent["betweenness_centrality"])
cent["eig_zscore"] = zscore(cent["eigenvector_centrality"])

# Combined hub score
cent["hub_score"] = (
    cent["deg_zscore"] +
    cent["bet_zscore"] +
    cent["eig_zscore"]
)

# Hub genes: combined z-score > 2 (top ~2.5%)
cent["is_hub"] = cent["hub_score"] > 2.0

# Separate TF hubs from target hubs
tf_list = set(pd.read_csv("tf_genes.txt", header=None)[0])
cent["is_tf"] = cent["gene"].isin(tf_list)

tf_hubs = cent[cent["is_hub"] & cent["is_tf"]].sort_values(
    "hub_score", ascending=False
)

print(f"Total hubs: {cent['is_hub'].sum()}")
print(f"TF hubs: {len(tf_hubs)}")
print(f"Top 10 hub TFs:\n{tf_hubs.head(10)[['gene', 'hub_score']]}")
```

## Hub-Target Expression Correlation Validation

```python
# Load expression matrix
expr = pd.read_csv("expression_matrix.tsv", sep="\t", index_col=0)

# Load edges
edges = pd.read_csv("grn_edges_genie3.tsv", sep="\t")

# For each hub TF, compute Spearman correlation with its targets
correlation_results = []
for tf in tf_hubs["gene"].head(20):
    tf_targets = edges[edges["TF"] == tf]["target"].tolist()
    tf_expr = expr.loc[tf] if tf in expr.index else None
    if tf_expr is None:
        continue

    corr_values = []
    for target in tf_targets:
        if target in expr.index:
            corr = np.corrcoef(tf_expr, expr.loc[target])[0, 1]
            if not np.isnan(corr):
                corr_values.append(corr)

    if corr_values:
        correlation_results.append({
            "TF": tf,
            "n_targets": len(corr_values),
            "mean_corr": np.mean(np.abs(corr_values)),
            "pct_positive": sum(1 for c in corr_values if c > 0) / len(corr_values) * 100
        })

corr_df = pd.DataFrame(correlation_results)
corr_df = corr_df.sort_values("mean_corr", ascending=False)
corr_df.to_csv("hub_tf_target_correlation.tsv", sep="\t", index=False)

print("Hub TF — Target correlation summary:")
print(corr_df.head(10))
```

## TF Family Enrichment in Hubs

```python
# Load TF family annotation
tf_family = pd.read_csv("tf_annotation.tsv", sep="\t")

# Merge family info with hub list
hub_with_family = tf_hubs.merge(
    tf_family[["gene_id", "family"]],
    left_on="gene", right_on="gene_id", how="left"
)

family_counts = hub_with_family["family"].value_counts()
print(f"TF families in hubs (top 10):\n{family_counts.head(10)}")

# Check: are plant-specific TF families enriched in hubs?
plant_families = ["MYB", "WRKY", "NAC", "bHLH", "AP2/ERF", "bZIP",
                   "GRAS", "MADS-box", "C2H2", "HD-ZIP", "TCP"]
for fam in plant_families:
    in_hub = family_counts.get(fam, 0)
    total_in_genome = len(tf_family[tf_family["family"] == fam])
    if total_in_genome > 0:
        pct_hub = in_hub / total_in_genome * 100
        print(f"  {fam}: {in_hub}/{total_in_genome} in hubs ({pct_hub:.1f}%)")
```

## Plant-Specific Notes

- **Expected hub TF families**: In most plant expression datasets, MYB, WRKY, NAC,
  and bHLH families consistently produce the highest number of hub TFs due to their
  large family sizes and central regulatory roles.

- **Hub gene stability**: Rank hubs by "robustness" — run GENIE3 on bootstrapped
  subsets (80% samples, 10 iterations). A TF that appears as a hub in >= 8/10
  bootstrap runs is a "stable hub".

- **Polyploid hub analysis**: In allopolyploids, check whether homeologous TF pairs
  both appear as hubs or only one copy. Single-copy hubs in polyploids suggest
  regulatory neofunctionalization or subfunctionalization.

- **Stress vs. control hubs**: If data includes multiple conditions, run hub analysis
  separately. Condition-specific hubs are the most interesting for functional validation.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Zero hub TFs identified | Network too sparse or threshold too high | Lower Z-score cutoff to 1.5; check network edge count |
| All hubs are non-TF target genes | Undirected network or TF list not applied | Verify GENIE3 was run with TF regulator constraint |
| Correlation with targets is near zero | Expression data not normalized | Check for batch effects; use Spearman not Pearson |
| Hub family enrichment all "unknown" | TF family annotation incomplete | Supplement with PFAM domain scan |
