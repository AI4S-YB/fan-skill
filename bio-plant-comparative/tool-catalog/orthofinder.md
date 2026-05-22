# OrthoFinder

**Goal:** Infer orthogroups, orthologs, and species tree from protein sequences
**Best for:** Multi-species comparative genomics — the gold standard

## Prerequisites
- OrthoFinder 2.5+ (Python)
- Protein sequences (.faa) for each species
- DIAMOND recommended for speed

## Basic Usage

```bash
mkdir proteomes && cp species*.faa proteomes/
orthofinder -f proteomes/ -t 16 -a 8
```

## Key Outputs
| File | Content |
|------|---------|
| Orthogroups.tsv | Gene family assignments across species |
| Species_Tree/SpeciesTree_rooted.txt | Inferred species tree |
| Orthologues/ | One-to-one, one-to-many orthologs |

## Plant-Specific Notes
- Plant genomes have many paralogs (due to WGD) — OrthoFinder handles this well
- Single-copy orthologs are useful for species tree construction
- Output orthogroups can feed into CAFE for gene family expansion/contraction
