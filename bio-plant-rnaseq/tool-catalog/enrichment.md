# Functional Enrichment Analysis

**Goal:** Identify over-represented pathways/GO terms among DEGs or co-expression modules
**Best for:** After DE or WGCNA — adds biological interpretation

## Prerequisites
- R 4.0+, clusterProfiler (Bioconductor)
- DEG list with gene IDs

## Plant Reactome ORA

```r
library(clusterProfiler)

# Download Plant Reactome GMT for your species
# https://plantreactome.gramene.org/download.html

gmt <- read.gmt("plant_reactome_rice.gmt")
enrich <- enricher(gene = deg_genes, TERM2GENE = gmt)
dotplot(enrich)
```

## GO Enrichment (AgriGO / clusterProfiler)

```r
library(clusterProfiler)
library(org.Osativa.eg.db)  # Rice example

enrich_go <- enrichGO(
  gene = deg_genes,
  OrgDb = org.Osativa.eg.db,
  ont = "BP",
  pAdjustMethod = "BH"
)
```

## Key Decisions

| Scenario | Best Tool |
|----------|----------|
| Rice, maize, Arabidopsis → DEG enrichment | Plant Reactome first, then GO |
| Non-model crop → DEG enrichment | Mercator annotation → GO via AgriGO |
| WGCNA module → enrichment | GO + KEGG (whole module as background) |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No genes mapped" | Gene ID format mismatch | Check ID type (RAP-DB vs MSU for rice) |
| "No enrichment found" | Too few DEGs or wrong background | Reduce stringency or use all expressed genes as background |
