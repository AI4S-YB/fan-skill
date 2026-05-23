# PTM Analysis (Motif-x / Phosphorylation)

**Goal:** Analyze post-translational modifications (PTMs), particularly phosphorylation, in plant proteomics
**Best for:** Phosphoproteomics, acetylomics, and other PTM-focused experiments in plants

## Prerequisites

- MaxQuant output (Phospho (STY)Sites.txt for phosphorylation)
- Motif-x (https://motif-x.med.harvard.edu/)
- Python 3.6+ with pandas, numpy
- Kinase-substrate database (PhosphoSitePlus, Plant PTM Viewer)

## Phosphorylation Site Analysis

### Loading and Filtering Phosphosites

```python
import pandas as pd
import numpy as np

# Load MaxQuant phosphosite output
phos = pd.read_csv("Phospho (STY)Sites.txt", sep="\t")

# Filter to class I sites (localization probability > 0.75)
phos_filtered = phos[
    (phos["Localization prob"] > 0.75) &
    (phos["Reverse"] != "+") &
    (phos["Potential contaminant"] != "+")
]

# Extract the 15-mer sequence window
phos_filtered["window"] = phos_filtered["Sequence window"]

# Separate by residue
phos_serine = phos_filtered[phos_filtered["Amino acid"] == "S"]
phos_threonine = phos_filtered[phos_filtered["Amino acid"] == "T"]
phos_tyrosine = phos_filtered[phos_filtered["Amino acid"] == "Y"]

print(f"Total class I sites: {len(phos_filtered)}")
print(f"  Serine (pS): {len(phos_serine)}")
print(f"  Threonine (pT): {len(phos_threonine)}")
print(f"  Tyrosine (pY): {len(phos_tyrosine)}")
print(f"  pS/pT/pY ratio: {len(phos_serine)/len(phos_filtered):.1%}/"
      f"{len(phos_threonine)/len(phos_filtered):.1%}/"
      f"{len(phos_tyrosine)/len(phos_filtered):.1%}")
```

## Motif-x Analysis

```python
# Prepare sequences for Motif-x
# Motif-x format: central residue in uppercase, flanking in lowercase
def prepare_motifx(sequences, center_pos=7):
    """Prepare 15-mer sequences for Motif-x"""
    motifs = []
    for seq in sequences:
        if len(seq) >= 15:
            motif = seq[:center_pos].lower() + seq[center_pos].upper() + seq[center_pos+1:].lower()
            motifs.append(motif)
    return motifs

# Write to file
with open("phospho_sites_motifx.txt", "w") as f:
    for m in prepare_motifx(phos_serine["window"].tolist()):
        f.write(m + "\n")
```

### Motif-x Web Submission

1. Go to https://motif-x.med.harvard.edu/
2. Upload the motif file
3. Set parameters:
   - Foreground: uploaded file
   - Background: species proteome (or UniProt species database)
   - Central character: S/T/Y
   - Width: 15
   - Occurrence: 20
   - Significance: 1e-6
4. Submit and download results

## Kinase Enrichment Analysis

```python
# Map phosphosites to known kinase motifs
# Based on kinase-substrate relationships

kinase_motifs = {
    "SnRK2": ["R-x-x-S"],       # ABA signaling (plant-specific)
    "CPK/CDPK": ["R/K-x-x-S"],  # Calcium-dependent (plant)
    "MPK3/6": ["P-x-S-P"],      # MAP kinase (stress response)
    "CK2": ["S/T-x-x-E/D"],     # Casein kinase 2
    "SnRK1": ["M/L-x-R-x-x-S"], # Energy sensing (plant)
}

def match_kinase_motif(sequence, motif):
    import re
    motif_pattern = motif.replace("x", "[A-Z]").replace("-", "")
    return bool(re.search(motif_pattern, sequence))

# Count sites matching each kinase
for kinase, motifs in kinase_motifs.items():
    count = 0
    for seq in phos_serine["window"]:
        for m in motifs:
            if match_kinase_motif(seq, m):
                count += 1
                break
    print(f"{kinase}: {count} sites")
```

## Differential Phosphorylation Analysis

```r
# Similar to protein-level DE, but analyze phosphosites
library(limma)

# Filter phosphosites with >= 3 valid values in one condition
phos_filt <- phos_mat[rowSums(!is.na(phos_mat)) >= 3, ]
phos_imp <- impute.knn(phos_filt, k = 10)$data

# limma analysis (same as protein-level)
fit <- lmFit(phos_imp, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
phos_results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

# Key difference from protein-level:
# |log2FC| >= 0.58 is standard for phospho (not 1.0)
phos_results$significant <- phos_results$adj.P.Val < 0.05 &
                             abs(phos_results$logFC) > 0.58
```

## Plant-Specific Notes

- Plant phosphosite stoichiometry is generally lower than mammalian (more transient phosphorylation)
- Typical yield: 3,000-8,000 phosphosites in plant leaf, 5,000-15,000 in roots
- Plant-specific phosphomotifs differ from mammalian consensus (e.g., plant SnRK2 motif)
- Multi-phosphorylated peptides: common in plants, especially in disordered regions
- For non-model plants: use PhosphoSitePlus plant data (Arabidopsis, rice, maize) for reference

## PTM Crosstalk Analysis

```python
# Check for co-occurring PTMs on the same protein
# Load multiple PTM datasets
phos = pd.read_csv("Phospho (STY)Sites.txt", sep="\t")
acet = pd.read_csv("Acetyl (K)Sites.txt", sep="\t")
ubi = pd.read_csv("Ubiquitin (K)Sites.txt", sep="\t")

# Merge on protein ID
merged = phos[["Protein", "Positions within proteins"]].merge(
    acet[["Protein", "Positions within proteins"]],
    on="Protein", suffixes=("_phos", "_acet"), how="inner"
)

print(f"Proteins with both phosphorylation and acetylation: {len(merged)}")
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Motif-x "no significant motifs" | Too few sites or too strict parameters | Lower occurrence to 10, relax significance to 1e-4 |
| High percentage of non-localized sites | Poor fragmentation | Increase MS/MS parameters, use phospho-enriched samples |
| No sites for non-model species | Lack of reference | Use BLAST to map to model species motif database |
| Very low phospho identification | Enrichment failed | Check TiO2/IMAC column quality, add lactic acid |
