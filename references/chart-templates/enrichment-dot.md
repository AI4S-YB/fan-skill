# Enrichment Dot Plot

**For**: GO/KEGG enrichment results from clusterProfiler, topGO, or similar tools

## Data Format Required

clusterProfiler output or equivalent:

| Column | Type | Description |
|--------|------|-------------|
| ID | character | GO term ID or KEGG pathway ID |
| Description | character | Human-readable term/pathway name |
| GeneRatio | character | "k/n" — number of DEGs in term / total DEGs |
| BgRatio | character/factor | Background gene ratio |
| p.adjust | numeric | FDR-adjusted p-value |
| Count | integer | Number of DEGs mapped to this term |
| Category | character/factor | GO category (BP/MF/CC) or KEGG |

## ggplot2 Code

```r
library(ggplot2)
library(dplyr)
library(stringr)
source("theme/theme_plant_scientific.R")

# Load enrichment results
# enrich <- read.table("go_enrichment.txt", header = TRUE, sep = "\t")

# Parse GeneRatio to numeric
enrich <- enrich %>%
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    Description = factor(Description, levels = Description[order(GeneRatio_num)]),
    neg_log10_padj = -log10(p.adjust)
  )

# Select top terms (top 10-20 per category)
top_terms <- enrich %>%
  group_by(Category) %>%
  slice_min(p.adjust, n = 10) %>%
  ungroup()

# ---- Dot Plot ----
ggplot(top_terms, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num))) +
  geom_point(aes(size = Count, color = neg_log10_padj)) +
  scale_color_viridis_c(
    name = expression(-log[10](adjusted ~ italic(P))),
    option = "D",
    direction = -1
  ) +
  scale_size_continuous(name = "Gene count", range = c(2, 8)) +
  scale_x_continuous(labels = scales::percent) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Gene Ratio",
    y = ""
  ) +
  theme_plant_scientific() +
  theme(
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(0.5, "lines")
  )
```

## Bar Chart Alternative

```r
# For fewer terms (<30) where reading exact enrichment is important
ggplot(top_terms, aes(x = Count, y = reorder(Description, Count))) +
  geom_col(aes(fill = neg_log10_padj), width = 0.7) +
  scale_fill_viridis_c(option = "D", direction = -1) +
  labs(
    x = "Gene Count",
    y = "",
    fill = expression(-log[10](adjusted ~ italic(P)))
  ) +
  theme_plant_scientific() +
  theme(aspect.ratio = 1.2)
```

## Enrichment Network (cnetplot)

```r
# For showing gene-term relationships
library(enrichplot)
cnetplot(
  enrich_result,
  showCategory = 10,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "all",
  cex_label_gene = 0.5,
  cex_label_category = 0.7
)
```

## Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `top_terms` per category | 10 | Increase for more comprehensive view; decrease for main figure |
| Point size range | c(2, 8) | Adjust for data range; keep minimum >= 2 |
| `direction` (viridis) | -1 | -1 = light for low p, dark for significant; 1 = reverse |
| Category | GO BP/MF/CC | Separate facets make comparisons across categories clearer |

## Plant-Specific Notes

- Plant-specific ontologies: Plant Ontology (PO), Plant Trait Ontology (TO), Gene Ontology slim for plants
- KEGG pathway maps: plants have species-specific pathways (photosynthesis, phenylpropanoid, etc.)
- For MapMan bins in plants: use bar chart grouped by MapMan category instead of dot plot
- If enrichment is very broad (many significant GO terms): use semantic similarity (rrvgo) to reduce redundancy before plotting
- For multi-species comparison: consider a comparison dot plot with species on x-axis groups within categories
