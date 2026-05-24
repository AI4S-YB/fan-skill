# GENIE3 — Tree-Based Network Inference

**Goal:** Infer directed gene regulatory networks from expression data using
random forest regression trees to predict each gene's expression from TF expression.

**Best for:** >= 15 samples, TF annotation available, directed edge inference.

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| nTrees / n_estimators | 1000 | Quick exploration: 500; publication final: 2000; >20,000 genes: use GRNBoost2 | Each tree captures a random subset of relationships; more trees stabilize edge weights at the cost of linear runtime increase |
| K (candidate regulators) | sqrt(n_genes) | TF-poor species (<100 known TFs): use "all"; TF-rich species (>3000 TFs): use sqrt(n_TFs) | K controls the random subset of regulators tested at each tree node; "all" TFs ensures complete search for poorly annotated plants |
| threshold (edge filtering) | 0.001 | Dense network (>100K edges): raise to 0.005; sparse network (<5K edges): lower to 0.0005 | Post-hoc threshold controls final network density; use network connectedness metrics to calibrate |
| TF list source | PlantTFDB / PlantRegMap | Non-model species: use OrthoFinder with Arabidopsis TFs; specific family analysis: custom list from InterPro/Pfam | TF annotation quality directly determines network completeness; many plant species lack curated TF lists |

## Prerequisites

- R 4.0+ or Python 3.8+
- R packages: GENIE3 (Bioconductor)
- Python packages: arboreto, pandas, numpy, scikit-learn (for GRNBoost2)
- Input: normalized expression matrix (genes x samples), TF gene list

## Basic Usage (R)

```r
library(GENIE3)

# Load expression matrix: rows = genes, columns = samples
exprMat <- as.matrix(read.table("expression_matrix.tsv",
    header = TRUE, row.names = 1, sep = "\t"))

# Load TF list (one gene ID per line)
tf_list <- scan("tf_genes.txt", what = "character")

# Run GENIE3
weightMat <- GENIE3(
    exprMat,
    regulators = tf_list,           # Only TFs as regulators
    nTrees = 1000,                  # Number of random forest trees
    nCores = 4,                     # Parallel cores
    seed = 42,
    verbose = TRUE
)

# Convert to edge list
linkList <- getLinkList(weightMat, threshold = 0.001)
write.table(linkList, "grn_edges_genie3.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE)
```

## Basic Usage (Python — GRNBoost2, faster alternative)

```python
import pandas as pd
from arboreto.algo import grnboost2

# Load expression matrix
expr = pd.read_csv("expression_matrix.tsv", sep="\t", index_col=0)

# Load TF list
with open("tf_genes.txt") as f:
    tf_list = [line.strip() for line in f]

# Run GRNBoost2
network = grnboost2(
    expression_data=expr,
    tf_names=tf_list,
    verbose=True
)

# Filter and save
network = network[network['importance'] > 0.001]
network.to_csv("grn_edges_genie3.tsv", sep="\t", index=False)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| nTrees / n_estimators | 1000 | Balance accuracy vs. runtime; 500 for quick run |
| regulators | TF gene list | MUST provide for directed network |
| threshold | 0.001 | Lower = more edges (dense); higher = more sparse |
| nCores | available - 1 | Parallelize tree building |

## Plant-Specific Notes

- **TF family enrichment**: After inference, check that top-ranked edges involve
  known plant TF families (MYB, WRKY, NAC, bHLH, AP2/ERF). Absence is a red flag.
- **Polyploid homeologs**: If two homeologous TFs have nearly identical target sets,
  flag them as potential regulatory redundancy. Distinguish expression divergence.
- **Tissue-specific networks**: Run GENIE3 per tissue; TFs active in specific tissues
  (e.g., MADS-box in floral tissue) should show tissue-specific network topology.
- **Drought/stress conditions**: Plant TF networks shift dramatically under stress
  (WRKY, NAC, DREB families drive rewiring). Run networks per condition.
- **Low expression filter**: Filter genes with mean TPM < 1 before GENIE3.
  Low expression introduces noise that inflates false-positive edges.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "regulators not in exprMat" | TF gene IDs don't match matrix | Harmonize gene IDs (check version/annotation) |
| Memory error | Matrix too large (>20K genes) | Pre-filter to top 5000-10000 variable genes |
| Runtime > 2h | Too many genes + trees | Use GRNBoost2 (Python), reduce genes, or set nTrees=500 |
| All edges near zero | Expression variance too low | Check if data was log-transformed; increase variance filter |
| No known TF families in top edges | Wrong annotation or non-plant | Verify TF database source; check species ID mapping |
