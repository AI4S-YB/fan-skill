# Plant Proteomics Specific Considerations

## Plant Protein Extraction Challenges

### Secondary Metabolites
Plant tissues are rich in compounds that interfere with protein extraction:
- **Phenolics**: Oxidize and cross-link proteins during extraction
- **Polysaccharides**: Increase viscosity, interfere with SDS-PAGE
- **Pigments** (chlorophyll, anthocyanins): Absorb at UV wavelengths, interfere with protein assays
- **Proteases**: Released during tissue disruption, degrade proteins

### Recommended Extraction Methods

| Tissue Type | Method | Rationale |
|-------------|--------|-----------|
| Green leaf | TCA/acetone precipitation | Removes chlorophyll, inhibits proteases |
| Seed/grain | Phenol extraction | Lipids and storage proteins |
| Root | TCA/acetone + PVPP | PVPP binds phenolics |
| Fruit | Phenol/SDS | High sugar and pigment content |
| Woody stem | SDS + sonication | Tough cell walls need mechanical disruption |

### Protein Extraction Buffer Recipe (General Plant)

```
TCA/Acetone Precipitation:
1. Grind tissue in liquid N2
2. Add 10% TCA in acetone (-20 C), vortex
3. Precipitate -20 C overnight
4. Wash 3x with cold acetone
5. Air dry, resuspend in UA buffer (8M urea, 0.1M Tris-HCl pH 8.5)
```

## RuBisCO Depletion

RuBisCO (Ribulose-1,5-bisphosphate carboxylase/oxygenase) constitutes 30-50% of total protein in green plant tissues:

### Depletion Methods

1. **Immunodepletion**: Anti-RuBisCO antibody columns (available for Arabidopsis, rice, maize)
2. **PEG precipitation**: 15% PEG 4000 selectively precipitates RuBisCO
3. **Heat precipitation**: 75 C for 10 min (RuBisCO is heat-labile)
4. **Ammonium sulfate fractionation**: 35-45% saturation

### When to Deplete

- **Always** for leaf proteomics if interested in non-photosynthetic proteins
- **Not needed** for non-green tissues (root, seed, tuber, fruit)
- **Not needed** if RuBisCO itself is the target

## Plant Protein Databases

### Model Species (Complete UniProt Reference Proteomes)

| Species | UniProt ID | Protein Count |
|---------|-----------|---------------|
| Arabidopsis thaliana | UP000006548 | ~27,000 |
| Oryza sativa (japonica) | UP000059680 | ~35,000 |
| Zea mays | UP000007305 | ~40,000 |
| Glycine max | UP000008827 | ~46,000 |
| Solanum lycopersicum | UP000004994 | ~34,000 |
| Triticum aestivum | UP000019116 | ~110,000 |
| Brachypodium distachyon | UP000008810 | ~26,000 |
| Medicago truncatula | UP000002051 | ~45,000 |

### Non-Model Species Strategy

1. **RNA-seq guided**: Translate assembled transcripts (6-frame) as custom database
2. **Orthology-based**: Use OrthoFinder results to map to model species
3. **BLAST-based**: Use protein sequences from closely related species (taxonomy level: genus/family)

## Post-Translational Modifications in Plants

### Phosphorylation
- Most studied plant PTM
- Plant-specific phosphosites often differ from mammalian consensus sequences
- Enrichment: TiO2 + lactic acid improves specificity for plant samples
- Typical: 2,000-5,000 phosphosites identified per experiment

### Acetylation
- Highly abundant in chloroplast proteins (especially photosynthetic enzymes)
- Histone acetylation regulates flowering time and stress response
- Enrichment: anti-acetyllysine antibody

### Ubiquitination
- Protein degradation via ubiquitin-proteasome system
- K48-linked: degradation signal
- K63-linked: signaling
- Enrichment: diGly remnant (K-GG) antibody after trypsin digestion

### Plant-Specific PTMs
- **Tyrosine nitration**: Nitric oxide signaling in stress response
- **Carbonylation**: Oxidative stress marker
- **SUMOylation**: Stress response, flowering time regulation

## Quantitative Strategies for Plant Proteomics

| Strategy | Pros | Cons | Best For |
|----------|------|------|----------|
| Label-free (LFQ) | Low cost, simple | Missing values, run-to-run variation | Large-scale screening, tissue comparisons |
| TMT (6-18 plex) | No missing values, high precision | Expensive, ratio compression | Precise quantification, time series |
| SILAC (in vivo) | Gold standard accuracy | Not feasible for plants (auxotrophic) | Not applicable for most plants |
| SILAC (in vitro) | Good accuracy | Only relative quantification | Spike-in standard |
| PRM/MRM | Highest specificity, absolute quant | Requires prior knowledge, low throughput | Biomarker validation |

## Plant-Specific Contaminants

### Common Background Proteins

| Protein | Source | Mitigation |
|---------|--------|------------|
| RuBisCO large subunit | Chloroplast | Depletion (see above) |
| RuBisCO small subunit | Chloroplast (nuclear-encoded) | Depletion |
| LHCB (Light Harvesting Complex) | Thylakoid membrane | Depletion or subcellular fractionation |
| Seed storage proteins | Seeds/grains | Selective extraction or depletion |
| PR proteins (PR1-PR17) | Pathogen-stressed tissue | Account for in experimental design |

### cRAP Database for Plants

Add these plant contaminants to the search database:
- Bovine serum albumin (BSA)
- Trypsin (porcine/bovine)
- Human keratins
- RuBisCO from non-target species (if using mixed samples)
