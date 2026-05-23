---
skill_name: bio-plant-chipseq
version: "1.0"
category: plant-bioinformatics
analysis_type: chip_sequencing
description: >
  Plant ChIP-seq data analysis pipeline.
  Covers peak calling for transcription factors and histone modifications,
  differential binding analysis, motif enrichment, peak annotation,
  quality control, and visualization for plant ChIP-seq data.
platform: Linux
language: [bash, R, python]
---

# bio-plant-chipseq

## Overview

This Skill provides a structured decision matrix and tool catalog for analyzing plant ChIP-seq data. It covers the full workflow from raw FASTQ preprocessing through peak calling, differential binding analysis, and functional annotation.

## When to Use This Skill

- You have plant ChIP-seq data and need to call peaks for transcription factors
- You have plant ChIP-seq data for histone modifications
- You need to perform differential binding analysis between conditions
- You need to annotate ChIP-seq peaks to nearby genes
- You need motif enrichment analysis of ChIP-seq peak regions
- You need quality control and visualization of ChIP-seq data

## Decision Matrix

The `decision-matrix.yaml` file encodes rules for selecting analysis strategies based on experimental conditions:

| Analysis Step          | Key Conditions                              | Recommendation                |
|------------------------|---------------------------------------------|-------------------------------|
| Peak Calling           | Transcription factor ChIP-seq               | MACS2 narrow peak             |
| Peak Calling           | Histone modification ChIP-seq               | MACS2 broad peak              |
| Differential Binding   | >= 2 biological replicates                  | DiffBind + DESeq2             |
| Peak Annotation        | Plant genome annotation                     | ChIPseeker + plant OrgDb      |

## Tool Catalog

| File                              | Description                                      |
|-----------------------------------|--------------------------------------------------|
| tool-catalog/preprocessing.md     | Quality control, read alignment (Bowtie2/BWA)    |
| tool-catalog/peak-calling.md      | MACS2 for TF and histone peaks                   |
| tool-catalog/diff-binding.md      | DiffBind, DESeq2 for differential binding        |
| tool-catalog/motif-enrichment.md  | MEME-ChIP, Homer for motif discovery            |
| tool-catalog/annotation.md        | ChIPseeker peak-to-gene annotation               |
| tool-catalog/visualization.md     | IGV tracks, heatmaps, profile plots             |

## References

- `references/chipseq-plant-special.md` — Plant-specific considerations for ChIP-seq
- `analyst-notebook.md` — Detailed Chinese-language analysis notebook

## Fallback

If no decision rule matches, the fallback rule delegates to expert review.
