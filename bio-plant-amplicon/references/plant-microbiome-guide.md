# Plant Microbiome Amplicon Guide

## Why Plant Microbiome Analysis Has Special Requirements

| Aspect | Human Microbiome | Plant Microbiome |
|--------|-----------------|------------------|
| Host DNA interference | Negligible | Major problem (chloroplast/mitochondria 16S) |
| Biomass | High (gut: 10^11 CFU/g) | Variable (soil: high; leaf: very low) |
| Spatial structure | Body sites distinct | Rhizosphere → endosphere gradient |
| Temporal dynamics | Stable in adults | Strongly seasonal + developmental stages |
| Environmental influence | Diet-driven | Soil type > host genotype |
| Reference databases | Gut-focused | Need soil/rhizosphere-specific references |

## Compartment-Specific Recommendations

### Bulk Soil
- Highest diversity and biomass
- Standard extraction and library prep
- Use soil-specific extraction kits (e.g., PowerSoil)
- Rarefaction depth: 10,000-30,000 reads/sample

### Rhizosphere
- Root-adhering soil (1-2mm from root surface)
- Enriched in copiotrophs (fast-growing bacteria)
- Extract within 24h of sampling or flash-freeze
- Standard depth: 5000-10,000 reads/sample

### Rhizoplane
- Root surface after washing off rhizosphere soil
- Sonication recommended for detachment
- Low biomass; PNA clamps strongly recommended
- Depth: 2000-5000 reads/sample

### Endosphere (Root Interior)
- Surface-sterilized root tissue
- Very low biomass; plant DNA > 90%
- PNA clamps ESSENTIAL
- Multiple negative controls required
- Depth: 1000-3000 reads/sample (often lower)

### Phyllosphere (Leaf Surface)
- Leaf washing or whole-leaf extraction
- Lowest biomass of all compartments
- UV exposure reduces bacterial load
- May need to pool leaves from same plant
- Depth: 500-2000 reads/sample

## PNA Clamps for Plant Studies

PNA (peptide nucleic acid) clamps block amplification of host organelle DNA:

| PNA Clamp | Target | Sequence |
|-----------|--------|----------|
| pPNA | Plant plastid 16S | GGCTCAACCCTGGACAG |
| mPNA | Plant mitochondria 18S | GGCAAGTCACCCTCCCA |

Typical PNA concentration: 0.5-1.0 uM in PCR reaction.
Add PNA to PCR master mix before template DNA.

## Negative Controls — Non-Negotiable

For plant microbiome studies, include these controls in EVERY run:

1. **DNA extraction blank**: Kit reagents processed without sample (1 per extraction batch)
2. **PCR negative**: PCR master mix without template (1 per PCR plate)
3. **Field blank**: Empty tube opened at field site during sampling (1 per sampling trip)

How to use controls:
- Sequenced controls at low depth (they should produce few reads)
- If a taxon is more abundant in samples than in controls, keep it
- If a taxon has similar abundance in samples and controls, flag it as potential contaminant
- Use `decontam` R package for statistical contaminant identification

## Plant-Specific Taxonomic Groups

### Commonly Enriched in Rhizosphere
- Proteobacteria (Alpha-, Gamma-)
- Actinobacteria
- Bacteroidetes
- Firmicutes (Bacillus, Paenibacillus)

### Common Plant Growth-Promoting Rhizobacteria (PGPR)
- Pseudomonas fluorescens group
- Bacillus subtilis group
- Azospirillum spp.
- Rhizobium spp. (legumes)
- Streptomyces spp.

### Common Fungal Endophytes
- Ascomycota (Pleosporales, Hypocreales)
- Basidiomycota (Sebacinales — including Serendipita/Piriformospora)

### Common Plant Pathogens Detected by ITS
- Fusarium oxysporum complex
- Alternaria spp.
- Rhizoctonia solani
- Phytophthora spp. (oomycete — not detected by ITS alone)

## Statistical Considerations

### Compositional Data
- All amplicon data is compositional (relative, not absolute)
- Use CLR (centered log-ratio) transform before parametric tests
- Avoid standard t-tests or ANOVA on relative abundance data
- ANCOM-BC and ALDEx2 are compositionally aware

### PCoA Analysis Interpretation
- PCoA axes rarely capture > 30% variance in microbiome data
- Do not over-interpret PCoA1 vs PCoA2 when combined variance < 20%
- Always report axis variance explained in figure captions
- Run PERMANOVA: report pseudo-F, R-squared, and p-value
- Run PERMDISP to check if significance is driven by location or dispersion
