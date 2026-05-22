---
name: bio-plant-population
description: >
  Plant population genetics — PCA, ADMIXTURE, Fst, phylogenetic trees.
  Characterize genetic structure, detect subpopulations, and quantify
  differentiation. Dual-mode decision system.
tool_type: mixed
primary_tool: plink2
workflow: true
---

# Plant Population Genetics

Genetic structure analysis for plant populations.

## Decision Modes

```yaml
decision_mode: hybrid   # rule | expert | hybrid
```

## Analysis Flow

### Step 1: LD Pruning
Prune linked SNPs before PCA/ADMIXTURE.
```bash
plink --bfile ${INPUT} --indep-pairwise 50 5 0.2 --out pruned
plink --bfile ${INPUT} --extract pruned.prune.in --make-bed --out input_ld
```

### Step 2: PCA
```bash
plink --bfile input_ld --pca 10 --out pca_result
```
→ Visualize PC1 vs PC2, check variance explained per PC.

### Step 3: ADMIXTURE (ancestry estimation)
Run K=1 to K=10, check CV error for best K.
```bash
for K in {1..10}; do
    admixture --cv input_ld.bed $K | tee admixture_${K}.log
done
```

### Step 4: Fst (pairwise differentiation)
If populations are known:
```bash
plink --bfile input_ld --fst --within population_assignments.txt --out fst_result
```

### Step 5: Phylogenetic Tree (NJ)
```bash
plink --bfile input_ld --distance square 1-ibs --out distance
Rscript -e 'library(ape); d=as.matrix(read.table("distance.mdist")); tr=nj(d); plot(tr, cex=0.5)'
```

### Step 6: Visualization
PCA scatter, ADMIXTURE bar plot, NJ tree, Fst heatmap.
See `tool-catalog/visualization.md`.

## QC Checkpoints

| Checkpoint | Metric | Expectation |
|------------|--------|-------------|
| After LD pruning | SNPs retained | > 50% of original |
| After PCA | PC1-3 cumulative | > 15% for structured populations |
| After ADMIXTURE | CV error | Converges (delta < 0.01 across last 3 K values) |

## References
- `bio-plant-infra/references/species-cheatsheet.md`
- `references/popgen-plant-special.md`
