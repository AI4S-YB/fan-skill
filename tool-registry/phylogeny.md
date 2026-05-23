# Phylogenetic Tree (Neighbor-Joining)

**Goal:** Visualize genetic relatedness among samples via distance-based tree
**Best for:** Understanding sample relationships complementary to PCA

## Prerequisites
- PLINK 1.9+
- R with `ape` package

## Basic Usage

```bash
# Calculate IBS distance matrix
plink --bfile input_ld --distance square 1-ibs --out distance

# Build NJ tree in R
Rscript -e '
library(ape)
d <- as.matrix(read.table("distance.mdist"))
tr <- nj(d)
pdf("outputs/figures/nj_tree.pdf", width=14, height=10)
plot(tr, cex=0.5, main="NJ Tree (1-IBS distance)")
dev.off()
'
```

## Plant-Specific Notes
- PLINK outputs lower-triangular distance matrix — must convert to full matrix
- Inbred lines cluster tightly; outcrossing populations show more branching
- If topology contradicts PCA, trust PCA for population-level inference
