# TF Database — Transcription Factor Annotation for Plants

**Goal:** Identify which genes in the expression matrix encode transcription factors,
using public databases (PlantTFDB) for model species or homology-based inference
for non-model species.

**Best for:** Any GRN analysis — TF annotation is the prerequisite for directed
network inference.

## Prerequisites

- Plant genome annotation (GFF3 + protein FASTA)
- Gene ID mapping between genome and expression matrix
- Python 3.8+: pandas, biopython, requests (for API)
- For homology inference: BLAST+ or DIAMOND, reference proteome (Arabidopsis)

## PlantTFDB — Model Species

Supported species: Arabidopsis thaliana, Oryza sativa, Zea mays, Glycine max,
Solanum lycopersicum, Medicago truncatula, Populus trichocarpa, and more.

```python
import pandas as pd
import requests

# PlantTFDB API: download TF list for species
# Base URL: http://planttfdb.gao-lab.org/

# Method 1: Download pre-computed TF list
# Visit http://planttfdb.gao-lab.org/download.php
# Select species and family classification

# Method 2: Parse downloaded TF list
def load_planttfdb(tf_file, species="arabidopsis_thaliana"):
    """Load PlantTFDB TF list and map to expression matrix gene IDs."""
    tf_df = pd.read_csv(tf_file, sep="\t", header=None,
                        names=["gene_id", "family"])
    print(f"Total TFs: {len(tf_df)}")
    print(f"TF families: {tf_df['family'].nunique()}")
    print(f"Top families:\n{tf_df['family'].value_counts().head(10)}")
    return set(tf_df["gene_id"])

# Usage
tf_set = load_planttfdb("Ath_TF_list.txt")
print(f"Found {len(tf_set)} TF genes")

# Save for GENIE3/SCENIC
with open("tf_genes.txt", "w") as f:
    for tf in sorted(tf_set):
        f.write(f"{tf}\n")
```

## Homology-Based TF Inference — Non-Model Species

```python
import subprocess
import pandas as pd

# Step 1: Run DIAMOND BLASTP against Arabidopsis proteome
subprocess.run([
    "diamond", "blastp",
    "--db", "arabidopsis_tf.dmnd",    # Pre-built: only Arabidopsis TF proteins
    "--query", "species_proteins.faa",
    "--out", "tf_blast.tsv",
    "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore",
    "--max-target-seqs", "1",
    "--evalue", "1e-5",
    "--threads", "4"
])

# Step 2: Filter high-confidence TF homologs
blast = pd.read_csv("tf_blast.tsv", sep="\t",
    names=["query", "target", "pident", "evalue", "bitscore"])
hits = blast[
    (blast["pident"] >= 30) &
    (blast["evalue"] < 1e-10)
]

# Step 3: Map back to expression matrix gene IDs
# (adapt ID mapping based on your species)
hits["tf_family"] = hits["target"].apply(
    lambda x: x.split("|")[1] if "|" in x else "unknown"
)

hits.to_csv("tf_annotation_homologs.tsv", sep="\t", index=False)
print(f"Inferred {len(hits['query'].unique())} TF genes from "
      f"{hits['tf_family'].nunique()} families")

tf_genes = set(hits["query"].unique())
with open("tf_genes.txt", "w") as f:
    for tf in sorted(tf_genes):
        f.write(f"{tf}\n")
```

## Plant-Specific Notes

- **Major plant TF families** (expect these to dominate):
  MYB, WRKY, NAC, bHLH, AP2/ERF (ERF/DREB), bZIP, GRAS, MADS-box,
  C2H2-ZnF, HD-ZIP, TCP, SBP, ARF, B3, Trihelix, GATA, LBD.

- **Family size varies by species**: Arabidopsis ~2000 TFs (~7% of genes),
  rice ~2400, maize ~3300 (due to genome expansion). Fewer than 500 TFs
  in your prediction likely indicates incomplete annotation.

- **Polyploid species**: Expect TF family count proportional to ploidy level.
  Homeologous TF pairs usually belong to the same family. Flag cases where
  only one homeolog is annotated as TF.

- **One-to-many orthology**: A non-model gene may BLAST to multiple Arabidopsis
  TFs. Accept all hits but flag for manual curation. In network inference,
  this creates "super-regulons" that may need splitting.

- **Orphan plant TFs**: Some plant-specific TFs have no clear Arabidopsis ortholog.
  These will be missed by homology inference. Check for InterPro domain presence
  (PFAM TF domains: PF00249 Myb_DNA-binding, PF03110 SBP, etc.) as complement.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| No TFs found | Gene IDs don't match DB naming | Check gene ID format (e.g., AT1G01030 vs AT1G01030.1) |
| Too few TFs (<100) | Missing annotations or wrong genome version | Try latest PlantTFDB or supplement with PFAM scan |
| DIAMOND no hits | Query sequences are nucleotide, not protein | Ensure input is protein FASTA (.faa), not CDS (.fna) |
| All TFs one family | BLAST threshold too permissive | Increase pident cutoff or filter by reciprocal best hit |
| TF list too large (>5000) | Including TF cofactors and chromatin remodelers | Filter to known DNA-binding domain families only |
