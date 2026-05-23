# MaxQuant (LFQ / TMT)

**Goal:** Protein identification and quantification from DDA (Data-Dependent Acquisition) MS data
**Best for:** Label-free quantification (LFQ) and TMT-based quantification in plant proteomics

## Prerequisites

- MaxQuant (https://www.maxquant.org/)
- .NET framework (Windows) or Mono (Linux)
- Raw MS files (.raw for Thermo, .mzML for other vendors)
- Plant protein FASTA database
- Contaminant database (included with MaxQuant)

## Basic Usage (LFQ)

```bash
# Linux command-line
maxquant mqpar.xml

# Typical configuration (mqpar.xml)
# - Raw files: all .raw files in experiment
# - FASTA: uniprot_arabidopsis.fasta + contaminants.fasta
# - Enzyme: Trypsin/P
# - Variable mods: Oxidation (M), Acetyl (Protein N-term)
# - Fixed mods: Carbamidomethyl (C)
# - LFQ: enabled, min ratio count = 2
# - Match between runs: enabled (0.7 min window)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| Enzyme | Trypsin/P | Standard for proteomics |
| Missed cleavages | 2 | Allows for incomplete digestion |
| Variable mods | Oxidation (M), Acetyl (N-term) | Common biological modifications |
| Fixed mods | Carbamidomethyl (C) | IAA alkylation |
| PSM FDR | 0.01 | 1% at PSM level |
| Protein FDR | 0.01 | 1% at protein level |
| LFQ min ratio count | 2 | At least 2 ratio counts for quantification |
| Match between runs | Enabled (0.7 min window) | Reduces missing values significantly |

## TMT Mode

```xml
<!-- mqpar.xml TMT section -->
<parameter name="IsobaricLabels">
  <string>TMTpro16plex</string>
</parameter>
<parameter name="Reporter ion MS2" />
<parameter name="Min. reporter PIF">0.75</parameter>
```

## Plant-Specific Notes

- Use species-specific FASTA from UniProt; for non-model plants, append RNA-seq translations
- Add plant-specific contaminants (RuBisCO, LHC, seed storage proteins) to contaminant list
- Match between runs is critical for plants — missing values are common due to secondary metabolite interference
- For polyploids: expect protein groups (not uniquely distinguishable proteins)
- Oxidation (M) is more common in plants due to higher oxidative stress
- Consider Phospho (STY) as variable mod if analyzing phosphoproteomics

## Output Files

| File | Content | Use |
|------|---------|-----|
| proteinGroups.txt | Protein-level quantification | Main results table |
| peptides.txt | Peptide-level quantification | Quality control |
| evidence.txt | MS/MS evidence | Detailed QC |
| Phospho (STY)Sites.txt | Phosphosite quantification | Phosphoproteomics (if enabled) |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No proteins identified" | Wrong FASTA or species mismatch | Verify database contains target species |
| "Out of memory" | Too many raw files | Run in batches of 10-20 files |
| Very high missing values | Match between runs disabled | Enable with 0.7 min window |
| "Mass recalibration failed" | Poor instrument calibration | Check if lock mass was enabled |
| Thousands of identifications lost | Too strict FDR | Relax protein FDR to 0.05 for exploratory |
