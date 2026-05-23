---
name: bio-plant-atacseq
version: "1.0"
category: plant-bioinformatics
analysis_type: atac_sequencing
description: >
  Plant ATAC-seq data analysis pipeline.
  Covers peak calling with MACS2, differential accessibility analysis
  with DiffBind+DESeq2, TF footprinting with TOBIAS, motif enrichment
  analysis, quality control, and visualization for plant ATAC-seq data.
platform: Linux
language: [bash, R, python]
---

# bio-plant-atacseq

## Overview

This Skill provides a structured decision matrix and tool catalog for analyzing plant ATAC-seq data. It covers the full workflow from raw FASTQ preprocessing through peak calling, differential accessibility analysis, transcription factor footprinting, and motif analysis.

## When to Use This Skill

- You have plant ATAC-seq data and need to identify open chromatin regions
- You need to perform peak calling for ATAC-seq using MACS2
- You need differential accessibility analysis between conditions
- You need transcription factor footprinting analysis
- You need motif enrichment analysis of accessible regions
- You need quality control and visualization of ATAC-seq data

## Decision Matrix

The `decision-matrix.yaml` file encodes rules for selecting analysis strategies based on experimental conditions:

| Analysis Step | Key Conditions | Recommendation |
|--------------|----------------|----------------|
| Peak Calling | ATAC-seq data | MACS2 --nomodel --shift -100 --extsize 200 |
| Differential Accessibility | >= 2 replicates | DiffBind + DESeq2 |
| Footprinting | TF analysis | TOBIAS footprinting |
| Motif Analysis | Peak summit sequences | MEME-ChIP / HOMER |

## Tool Catalog

| File | Description |
|------|-------------|
| tool-catalog/preprocessing.md | Quality control, read trimming, alignment (Bowtie2/BWA) |
| tool-catalog/peak-calling.md | MACS2 ATAC-seq mode peak calling |
| tool-catalog/diff-accessibility.md | DiffBind, DESeq2 for differential accessibility |
| tool-catalog/footprinting.md | TOBIAS TF footprinting analysis |
| tool-catalog/motif-analysis.md | MEME-ChIP, HOMER motif discovery |
| tool-catalog/visualization.md | IGV tracks, heatmaps, accessibility profiles |

## References

- `references/atacseq-plant-special.md` -- Plant-specific considerations for ATAC-seq
- `analyst-notebook.md` -- Detailed Chinese-language analysis notebook

## Fallback

If no decision rule matches, the fallback rule delegates to expert review.
