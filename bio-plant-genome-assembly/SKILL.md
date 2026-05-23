---
name: bio-plant-genome-assembly
version: "1.0"
category: plant-bioinformatics
analysis_type: genome_assembly
description: >
  Plant genome assembly pipeline.
  Covers long-read assembly with hifiasm (PacBio HiFi) and Flye (ONT),
  hybrid assembly, polishing with Medaka/gcpp, Hi-C scaffolding,
  quality assessment with BUSCO and Merqury, and visualization.
platform: Linux
language: [bash, R, python]
---

# bio-plant-genome-assembly

## Overview

This Skill provides a structured decision matrix and tool catalog for assembling plant genomes. It covers the full workflow from sequencing data through assembly, polishing, scaffolding, and quality assessment, with specific guidance for the unique challenges of plant genomes (polyploidy, repeat content, large genome size).

## When to Use This Skill

- You have PacBio HiFi or ONT long-read data for genome assembly
- You need a de novo plant genome assembly
- You have both long and short reads for hybrid assembly
- You need to polish a draft assembly
- You have Hi-C data for chromosome-level scaffolding
- You need to assess assembly quality with BUSCO and Merqury

## Decision Matrix

The `decision-matrix.yaml` file encodes rules for selecting assembly strategies based on sequencing data type:

| Analysis Step | Key Conditions | Recommendation |
|--------------|----------------|----------------|
| Assembly | PacBio HiFi | hifiasm |
| Assembly | ONT | Flye |
| Assembly | Long + Short reads | Hybrid assembly |
| Polishing | ONT reads | Medaka |
| Polishing | PacBio HiFi | gcpp |
| QC | Any assembly | BUSCO + Merqury |
| Scaffolding | Hi-C data available | YAHS or SALSA Hi-C |

## Tool Catalog

| File | Description |
|------|-------------|
| tool-catalog/hifiasm.md | hifiasm for PacBio HiFi assembly |
| tool-catalog/flye.md | Flye for ONT assembly |
| tool-catalog/polishing.md | Medaka, gcpp for assembly polishing |
| tool-catalog/qc.md | BUSCO completeness + Merqury k-mer evaluation |
| tool-catalog/scaffolding.md | YAHS, SALSA Hi-C scaffolding |
| tool-catalog/visualization.md | Assembly visualization, Bandage plots |

## References

- `references/assembly-plant-special.md` -- Plant-specific assembly considerations
- `analyst-notebook.md` -- Detailed Chinese-language analysis notebook

## Fallback

If no decision rule matches, the fallback rule delegates to expert review.
