# KASP Marker Design

**Goal:** Design Kompetitive Allele-Specific PCR assays for SNP genotyping
**Best for:** Low-to-medium throughput SNP genotyping in breeding

## Prerequisites
- SNP position and flanking sequence (>=100 bp each side)
- Primer3 for primer design

## Design Process

1. Extract flanking sequence (+-200 bp around SNP)
2. Design two allele-specific forward primers (SNP at 3' end) + common reverse
3. Add FRET cassette tails (FAM and HEX) to forward primers
4. Check specificity against reference genome

## Plant-Specific Notes
- Polyploid species: the SNP must distinguish between subgenomes
- KASP works well for 1-100 markers in breeding populations
- For >100 markers, consider targeted GBS or SNP chip instead
