---
name: bio-plant-methylation
version: "1.0"
category: plant-bioinformatics
analysis_type: dna_methylation
description: >
  Plant DNA methylation data analysis pipeline.
  Covers WGBS and RRBS alignment with Bismark, methylation calling,
  DMR detection with DSS, methylation-expression correlation analysis,
  quality control, and visualization for plant methylome data.
platform: Linux
language: [bash, R, python]
---

# bio-plant-methylation

## Overview

This Skill provides a structured decision matrix and tool catalog for analyzing plant DNA methylation data. It covers the full workflow from raw FASTQ preprocessing through alignment, methylation calling, differentially methylated region (DMR) detection, and correlation analysis with gene expression.

## When to Use This Skill

- You have plant WGBS or RRBS data and need methylation analysis
- You need Bismark alignment and methylation extraction
- You need DMR calling between conditions
- You need to integrate methylation data with gene expression data
- You need methylation data quality control and visualization
- You are working with plant-specific methylation contexts (CG, CHG, CHH)

## Decision Matrix

The `decision-matrix.yaml` file encodes rules for selecting analysis strategies based on experimental protocol:

| Analysis Step | Key Conditions | Recommendation |
|--------------|----------------|----------------|
| Alignment | WGBS protocol | Bismark WGBS alignment + methylation extraction |
| Alignment | RRBS protocol | Bismark RRBS mode |
| DMR Calling | >= 2 replicates | DSS DMR detection |
| Correlation | Has expression data | Methylation-expression correlation |
| Context Analysis | Plant genomes | CG, CHG, CHH separation |

## Tool Catalog

| File | Description |
|------|-------------|
| tool-catalog/preprocessing.md | Quality control, read trimming for bisulfite data |
| tool-catalog/bismark.md | Bismark alignment and methylation extraction |
| tool-catalog/methylation-calling.md | Methylation level calling per cytosine |
| tool-catalog/dmr-calling.md | DSS for differentially methylated regions |
| tool-catalog/correlation-analysis.md | Methylation-gene expression correlation |
| tool-catalog/visualization.md | Methylation level heatmaps, browser tracks |

## References

- `references/methylation-plant-special.md` -- Plant-specific methylation considerations
- `analyst-notebook.md` -- Detailed Chinese-language analysis notebook

## Fallback

If no decision rule matches, the fallback rule delegates to expert review.
