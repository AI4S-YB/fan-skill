---
skill_name: bio-plant-small-rna
version: "1.0"
category: plant-bioinformatics
analysis_type: small_rna_sequencing
description: >
  Plant small RNA sequencing data analysis pipeline.
  Covers miRNA prediction, target gene prediction, degradome analysis,
  differential expression, and quality control for plant small RNA-seq data.
platform: Linux
language: [bash, R, python]
---

# bio-plant-small-rna

## Overview

This Skill provides a structured decision matrix and tool catalog for analyzing plant small RNA sequencing data. It covers the full workflow from raw FASTQ preprocessing through miRNA identification, target prediction, and differential expression analysis.

## When to Use This Skill

- You have plant small RNA-seq data (FASTQ) and need to identify known/novel miRNAs
- You need to predict miRNA target genes in plants
- You need to perform differential expression analysis of small RNAs
- You need to analyze degradome sequencing data for miRNA target validation
- You need quality control and visualization of small RNA-seq data

## Decision Matrix

The `decision-matrix.yaml` file encodes rules for selecting analysis strategies based on experimental conditions:

| Analysis Step          | Key Conditions                              | Recommendation                    |
|------------------------|---------------------------------------------|-----------------------------------|
| miRNA Prediction       | Model species detected                      | miRDeep2 + miRBase known miRNAs   |
| miRNA Prediction       | Non-model / unknown species                 | miRDeep2 novel discovery          |
| Target Prediction      | Plant miRNA targets                         | psRNATarget                       |
| Differential Expression| >= 3 biological replicates                  | DESeq2 on miRNA counts            |

## Tool Catalog

| File                              | Description                                      |
|-----------------------------------|--------------------------------------------------|
| tool-catalog/preprocessing.md     | Quality control, adapter trimming, length filter |
| tool-catalog/mirna-prediction.md  | miRDeep2, miRBase integration                    |
| tool-catalog/target-prediction.md | psRNATarget, TargetFinder                        |
| tool-catalog/degradome.md         | CleaveLand, degradome-seq analysis              |
| tool-catalog/differential-mirna.md| DESeq2, edgeR for miRNA counts                  |
| tool-catalog/visualization.md     | Volcano plots, heatmaps, miRNA structure        |

## References

- `references/small-rna-plant-special.md` — Plant-specific considerations for small RNA analysis
- `analyst-notebook.md` — Detailed Chinese-language analysis notebook with step-by-step protocols

## Fallback

If no decision rule matches, the fallback rule delegates to expert review.
