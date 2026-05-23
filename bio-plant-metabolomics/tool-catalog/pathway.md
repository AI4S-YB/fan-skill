# Pathway Mapping

**Goal:** Map identified metabolites to plant metabolic pathways and perform enrichment analysis
**Best for:** Interpreting metabolomics results in a biological context

## Prerequisites
- Metabolite identifiers or names (from annotation)
- PlantCyc and/or KEGG database access
- R 4.0+ or Python 3 with relevant packages

## PlantCyc Pathway Mapping

### Access PlantCyc

```bash
# PlantCyc is available via:
# 1. Web: https://plantcyc.org/
# 2. API: https://websvc.plantcyc.org/
# 3. Download: Pathway Tools + plantcyc database
```

### Map Metabolites to Pathways (Python)

```python
import requests
import pandas as pd

# Load identified metabolites
metabolites = pd.read_csv("identified_metabolites.csv")

# Query PlantCyc API for pathway assignments
pathway_results = []
for _, row in metabolites.iterrows():
    compound_name = row['compound_name']
    # PlantCyc API: search compound
    resp = requests.get(
        f"https://websvc.plantcyc.org/apixml?fn=search&id=PLANT&query={compound_name}"
    )
    # Parse and store pathway info
    if resp.status_code == 200:
        pathway_results.append({
            'compound': compound_name,
            'pathways': parse_pathways(resp.text)
        })

# Create pathway-metabolite matrix for enrichment
pathway_df = pd.DataFrame(pathway_results)
```

### KEGG Plant Pathway Mapping

```r
library(KEGGREST)

# Get list of plant pathways
plant_pathways <- keggList("pathway", "ko")
plant_paths <- grep("Metabolic|Biosynthesis", plant_pathways, value = TRUE)

# Map compounds to pathways
for (cpd in identified_compounds) {
    # Search KEGG compound
    result <- keggFind("compound", cpd)
    pathways <- keggLink("pathway", names(result))
    # Store mapping
}
```

## Pathview: KEGG Pathway Visualization

```r
library(pathview)

# Visualize metabolite changes on KEGG pathway
# Requires: KEGG pathway ID + named vector of log2 fold changes

# Prepare data
logFC <- significant_metabolites$logFC
names(logFC) <- significant_metabolites$kegg_id

# Generate pathway view
pathview(
  cpd.data = logFC,
  pathway.id = "00941",        # Flavonoid biosynthesis
  species = "ko",              # KEGG Orthology (cross-species)
  kegg.dir = "kegg_pathways",
  low = list(gene = "blue", cpd = "blue"),
  mid = list(gene = "grey", cpd = "grey"),
  high = list(gene = "red", cpd = "red"),
  bins = list(gene = 10, cpd = 10),
  limit = list(gene = 2, cpd = 2)
)

# For plant-specific: use "ath" (Arabidopsis thaliana) or "osa" (Oryza sativa)
pathview(cpd.data = logFC, pathway.id = "00941", species = "ath")
```

## MetaboAnalystR (Comprehensive Web Alternative)

```r
# MetaboAnalystR provides web-based pathway analysis
# Download and install:
# install.packages("MetaboAnalystR", repos = "https://www.metaboanalyst.ca/MetaboAnalystR")

library(MetaboAnalystR)

# Create mSet object and run pathway analysis
mSet <- InitDataObjects("conc", "pathora", FALSE)
mSet <- Read.TextData(mSet, "metabolite_concentrations.csv", "rowu", "disc")
mSet <- CrossReferencing(mSet, "name")
mSet <- CreateMappingResultTable(mSet)
mSet <- SetKEGG.PathLib(mSet, "ath")  # Arabidopsis pathways
mSet <- CalculateOraScore(mSet, "hyperg", "pathway")
```

## Key Plant Metabolic Pathways

| Pathway | KEGG ID | PlantCyc ID | Example Compounds |
|---------|---------|-------------|-------------------|
| Flavonoid biosynthesis | map00941 | PWY-6797 | Quercetin, kaempferol, anthocyanins |
| Phenylpropanoid biosynthesis | map00940 | PHENYL-PWY | Caffeic acid, coumaric acid, lignin precursors |
| Terpenoid backbone biosynthesis | map00900 | NONMEVIPP-PWY | IPP, DMAPP |
| Glucosinolate biosynthesis | map00966 | PWY-5267 | Glucoraphanin (Brassica-specific) |
| Carotenoid biosynthesis | map00906 | CAROTENOID-PWY | Beta-carotene, lycopene, lutein |
| Alkaloid biosynthesis | map00960 | PWY-5858 | Nicotine, tropane alkaloids (Solanaceae) |
| Starch and sucrose metabolism | map00500 | PWY-622 | Glucose, fructose, sucrose, starch |

## PlantCyc Enrichment Analysis

```r
# Hypergeometric test for PlantCyc pathway enrichment
# Input: list of pathway IDs for identified compounds + background

enrich_plantcyc <- function(hit_pathways, background_pathways, plantcyc_pathways) {
  results <- data.frame()
  total_compounds <- length(background_pathways)
  hit_compounds <- length(hit_pathways)

  for (pathway in names(plantcyc_pathways)) {
    pathway_size <- length(plantcyc_pathways[[pathway]])
    overlap <- length(intersect(hit_pathways, plantcyc_pathways[[pathway]]))

    if (overlap >= 2) {
      p_val <- phyper(overlap - 1, pathway_size,
                      total_compounds - pathway_size,
                      hit_compounds, lower.tail = FALSE)
      results <- rbind(results, data.frame(
        pathway = pathway,
        overlap = overlap,
        pathway_size = pathway_size,
        p_value = p_val
      ))
    }
  }
  results$q_value <- p.adjust(results$p_value, "BH")
  return(results[order(results$q_value), ])
}
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| species | KEGG organism code (ath, osa, zma, sly, etc.) |
| limit | Log fold change range for color mapping |
| bins | Number of color bins |

## Plant-Specific Considerations

- PlantCyc has 600+ pathways vs KEGG's ~150 plant metabolic pathways
- KEGG is biased toward model organisms; non-model crop species may have divergent pathways
- Many plant secondary metabolite pathways are species- or family-specific
- Combine multiple KEGG maps for related pathways (e.g., map00940+map00941+map00944 for phenylpropanoid network)
- For GMO or edited plants, use the closest model species pathway as reference

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No pathway found" | Metabolite not in PlantCyc/KEGG | Use CANOPUS chemical class instead |
| "Download error" | KEGG API rate limit | Add Sys.sleep(1) between queries |
| Few metabolites mapped | Most features unidentified | Only report pathway analysis on identified subset; note limitation |
| Mismatched pathways | Species-specific pathway names | Verify with literature for your species |
