# DIA-NN

**Goal:** Protein identification and quantification from DIA (Data-Independent Acquisition) MS data
**Best for:** DIA-based plant proteomics with deep neural network-based analysis

## Prerequisites

- DIA-NN (https://github.com/vdemichev/DiaNN)
- DIA raw files (.raw or .mzML) or converted files (.d)
- Spectral library (optional — library-free mode available)
- Plant protein FASTA database
- GPU recommended but not required

## Basic Usage

```bash
# Command line (library-free mode)
DiaNN \
  --f "raw/*.raw" \
  --lib "" \
  --fasta "arabidopsis.fasta" \
  --threads 32 \
  --out "diann_output.tsv" \
  --qvalue 0.01 \
  --matrices

# With spectral library
DiaNN \
  --f "raw/*.raw" \
  --lib "arabidopsis_speclib.tsv" \
  --fasta "arabidopsis.fasta" \
  --threads 32 \
  --out "diann_output.tsv" \
  --qvalue 0.01

# Generate report matrix
DiaNN --matrices --gen-report
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `--qvalue` | 0.01 | FDR threshold (1%) |
| `--threads` | 16-32 | Thread count for parallel processing |
| `--matrices` | Enabled | Output protein x sample matrix |
| `--fasta` | Species-specific | Protein sequence database |
| `--lib` | "" (empty) | Library-free mode (deep learning predicts spectra) |
| `--mass-acc` | 10 ppm (Orbitrap), 20 ppm (TOF) | Mass accuracy |
| `--missed-cleavages` | 1 | Trypsin missed cleavages |
| `--var-mods` | 1 (Oxidation M) | Variable modifications count |

## Plant-Specific Notes

- DIA-NN's library-free mode works well for well-annotated plant species (Arabidopsis, rice, maize)
- For non-model plants: build a spectral library from DDA data or use predicted spectra
- The deep learning predictor generates in silico spectra — works best for tryptic peptides
- Plant-specific modifications (phosphorylation on STY) require separate DIA-NN runs with appropriate mod settings
- Output matrices have very low missing value rate (< 5% typically) — major advantage for plant samples

## Spectral Library Generation

```bash
# If you have DDA data from the same species
# 1. Generate spectral library with FragPipe or Spectronaut
# 2. Convert to DIA-NN format
# 3. Use --lib to specify the spectral library

# Alternatively, use Prosit (deep learning prediction)
# to generate predicted spectral library
```

## Output Processing

```python
import pandas as pd

# Load DIA-NN output
data = pd.read_csv("diann_output.tsv", sep="\t")

# Filter to proteotypic genes
data_filt = data[data["Proteotypic"] == 1]

# Extract protein-level matrix
pg_matrix = data_filt.pivot_table(
    index="Genes", columns="Run", values="PG.MaxLFQ",
    aggfunc="max"
)

# Log2 transform
import numpy as np
pg_log2 = np.log2(pg_matrix + 1)
pg_log2.to_csv("protein_lfq_matrix.csv")
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No identifications" | Library-free mode with non-plant FASTA | Check FASTA format (no spaces in headers) |
| Memory error | Too many raw files simultaneously | Process files in batches |
| "FASTA too large" | Plant genome too large | Filter to canonical isoforms, or split chromosomes |
| Low identification rate | Wrong mass accuracy | Check MS instrument specs and adjust `--mass-acc` |
| GPU out of memory | Batch size too large for GPU | Reduce batch size or use CPU mode |
