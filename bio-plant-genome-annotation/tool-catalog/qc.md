# Annotation Quality Control (BUSCO)

**Goal:** Assess completeness of genome annotation using conserved single-copy orthologs
**Best for:** Post-annotation quality assessment of any plant genome

## Prerequisites
- BUSCO 5.0+
- Lineage dataset (e.g., embryophyta_odb10, eudicots_odb10, poales_odb10)
- Predicted protein sequences (FASTA)

## Download Lineage Dataset

```bash
# List available plant lineages
busco --list-datasets | grep -i "plant\|embryophyta\|viridiplantae\|eudicots\|poales"

# Download appropriate dataset
busco --download embryophyta_odb10
```

## Run BUSCO in Protein Mode

```bash
busco \
  -i predicted_proteins.fasta \
  -l embryophyta_odb10 \
  -o busco_annotation \
  -m proteins \
  --cpu 16 \
  --long
```

## Alternative Modes

```bash
# For nucleotide gene sequences (CDS)
busco -i predicted_cds.fasta -l embryophyta_odb10 -o busco_cds -m transcriptome --cpu 16

# For genome assembly (before annotation)
busco -i genome.fasta -l embryophyta_odb10 -o busco_genome -m genome --cpu 16
```

## Interpreting BUSCO Results

The summary file (`short_summary.txt`) contains:

```
C:90.5%[S:87.2%,D:3.3%],F:5.1%,M:4.4%,n:1614
```

| Symbol | Meaning | What to Check |
|--------|---------|---------------|
| C | Complete | Should be > 80% |
| S | Complete + single-copy | Ideally > 80% |
| D | Complete + duplicated | < 10% (higher = assembly redundancy or polyploid) |
| F | Fragmented | < 10% (higher = incomplete gene models) |
| M | Missing | < 10% (higher = annotation missed conserved genes) |
| n | Total BUSCO genes | Reference dataset size |

## Plant-Specific Lineage Datasets

| Dataset | Taxa Coverage | Use Case |
|---------|---------------|----------|
| embryophyta_odb10 | All land plants (n=1614) | Default for any plant annotation |
| viridiplantae_odb10 | Green plants including algae (n=425) | Algae + land plant comparisons |
| eudicots_odb10 | Eudicots (n=2326) | Arabidopsis, soybean, tomato, cotton |
| poales_odb10 | Grasses and allies (n=4896) | Rice, maize, wheat, sorghum, barley |
| solanales_odb10 | Nightshade family (n=5950) | Tomato, potato, pepper |
| fabales_odb10 | Legumes (n=5366) | Soybean, common bean, pea |
| brassicales_odb10 | Mustard family (n=5387) | Arabidopsis, Brassica species |

## Quality Thresholds

| Grade | Complete | Fragmented | Action |
|-------|----------|------------|--------|
| **Excellent** | > 95% | < 3% | Annotation ready for publication |
| **Good** | 80–95% | < 10% | Acceptable; note caveats in methods |
| **Fair** | 70–80% | < 15% | Consider adding more RNA-seq evidence |
| **Poor** | < 70% | > 15% | Re-do gene prediction with better evidence |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Dataset not found" | Lineage not downloaded | Run `busco --download <dataset>` |
| "Input too short" | Proteins contain stop codons | Use protein sequences, not CDS |
| "All genes missing" | Wrong lineage or non-plant proteins | Verify input is plant proteins |
| "CPU limit exceeded" | Too many threads requested | Reduce --cpu; BUSCO spawns per-dataset processes |
| High duplication (>20%) | Polyploid or assembly redundancy | Expected for polyploids; note in interpretation |
