# CheckM (MAG Quality Assessment)

**Goal:** Assess metagenome-assembled genome (MAG) completeness and contamination
**Best for:** QC filtering of MAGs before downstream analysis, selecting high-quality genomes

## Prerequisites

- CheckM (https://github.com/Ecogenomics/CheckM)
- MAG FASTA files (one per bin)
- Lineage-specific marker gene set (downloaded automatically)

## Basic Usage

```bash
# Lineage-specific workflow (more accurate but slower)
checkm lineage_wf \
  -t 16 \
  -x fa \
  bins/ \
  checkm_output/

# Taxonomy assignment
checkm taxonomy_wf \
  domain Bacteria \
  -t 16 \
  -x fa \
  bins/ \
  checkm_taxonomy/

# Quick assessment (faster but less accurate)
checkm tetra bins/ checkm_tetra.tsv
checkm coverage bins/ checkm_coverage.tsv

# Quality summary
checkm qa \
  checkm_output/lineage.ms \
  checkm_output/ \
  -o 2 \
  -t 16 \
  --tab_table \
  -f mag_quality_report.tsv
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `lineage_wf` | Default workflow | Most accurate completeness/contamination |
| `-t` | 16-32 | Thread count |
| `-x` | fa/fasta | MAG file extension |
| `-o 2` (qa) | Output all fields | Get full QC report |
| `reduced_tree` | Use if runtime is too long | Faster but slightly less accurate |

## MAG Quality Classification

Based on MIMAG standards:

| Quality | Completeness | Contamination | 16S/18S/ITS | tRNA count |
|---------|-------------|---------------|--------------|------------|
| **High** | >90% | <5% | Present (full) | >=18 |
| **Medium** | >=50% | <10% | Present (partial) | <18 |
| **Low** | <50% | - | Absent | - |

## Plant-Specific Notes

- Plant-associated bacteria often have reduced genomes (endosymbionts/intercellular) — 80% completeness may be the best achievable
- Check for plant organelle contamination (chloroplast/mitochondria) in MAGs
- Rhizobia and other plant symbionts may have large plasmids affecting completeness estimates
- Use GTDB-Tk after CheckM for taxonomic classification

## Filtering High-Quality MAGs

```bash
# Extract high-quality MAGs
awk -F'\t' 'NR>1 && $6>90 && $7<5 {print $1}' \
  mag_quality_report.tsv > high_quality_mags.txt

# Extract medium+ quality MAGs
awk -F'\t' 'NR>1 && $6>=50 && $7<10 {print $1}' \
  mag_quality_report.tsv > medium_quality_mags.txt

# Copy quality MAGs to new directory
while read mag; do
  cp bins/${mag}.fa quality_bins/
done < high_quality_mags.txt
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No marker genes found" | Contig too short, or archaeal MAG in Bacteria mode | Check domain, increase contig length cutoff |
| "Completeness > 100%" | Multiple copies of single-copy markers = contamination | Contamination is already >5%, exclude MAG |
| "Lineage workflow too slow" | Too many MAGs | Use `reduced_tree` flag or pre-filter by size |
| Strain heterogeneity warning | Multiple strains binned together | Report strain heterogeneity; may inflate completeness |
