# Plant-Specific eQTL Analysis Reference

Guidance specific to expression quantitative trait locus mapping in plant/crop species.

## Polyploid eQTL Mapping Complexity

Polyploidy is the defining challenge of plant eQTL studies. Unlike diploid species where
each gene has one genomic position, polyploids present a multifaceted problem:

### Subgenome Complexity

- **Allopolyploids** (wheat AABBDD, cotton AADD, Brassica napus AACC): Homeologous genes
  from different subgenomes may share high sequence similarity. A SNP in one subgenome
  can appear as a "trans-eQTL" for homeologs in other subgenomes.
- **Autopolyploids** (potato AAAA, alfalfa AAAA): Allele dosage rather than biallelic
  genotypes. Standard eQTL tools assume diploid genotypes and will mis-estimate effect sizes.

### Recommended Strategies

1. **Subgenome-specific mapping**: Assign SNPs and genes to subgenomes before analysis.
   Use subgenome-specific reference (e.g., IWGSC RefSeq v2.1 for wheat with A/B/D
   annotations). Run eQTL separately per subgenome.

2. **Homeolog-aware expression quantification**: Use tools like HomeoRoq or Salmon with
   subgenome-specific transcriptomes to get accurate per-homeolog expression.

3. **Cross-homeolog SNP filtering**: Remove SNPs where reads map to multiple subgenomes
   (MAPQ < 30 or alternate mapping score). These are the main source of false trans-eQTLs.

4. **Tetraploid considerations**: For autotetraploids, use tools that handle dosage
   (e.g., GWASpoly adaptation for eQTL) or use presence/absence variant coding.

### Species-Specific Polyploid Notes

| Species | Ploidy | Genome size | Subgenomes | eQTL recommendation |
|---------|--------|-------------|------------|---------------------|
| Wheat (bread) | Allohexaploid | ~17 Gb | A, B, D | Per-subgenome cis only; trans is unreliable |
| Wheat (durum) | Allotetraploid | ~12 Gb | A, B | Per-subgenome, moderate power |
| Cotton (upland) | Allotetraploid | ~2.5 Gb | At, Dt | Per-subgenome, good power |
| Brassica napus | Allotetraploid | ~1.1 Gb | A, C | Per-subgenome, good power |
| Potato | Autotetraploid | ~0.84 Gb | - | Use dosage-aware methods |
| Sugarcane | Auto-octoploid | ~10 Gb | - | Extremely complex; cis-only |
| Strawberry | Allooctoploid | ~0.8 Gb | A, B1, B2, C | Per-subgenome, small genomes |

## Tissue Atlas Availability for Crops

Publicly available expression atlases for interpreting tissue-specific eQTL:

### Major Crop Expression Atlases

| Species | Resource | URL | Tissues | Notes |
|---------|----------|-----|---------|-------|
| Rice (Oryza sativa) | Rice eFP Browser | bar.utoronto.ca/efprice | 50+ | Best plant atlas; developmental series |
| Rice | RiceXPro | ricexpro.dna.affrc.go.jp | 30+ | Field-grown samples; stress treatments |
| Maize (Zea mays) | Maize Atlas / qTeller | qteller.maizegdb.org | 30+ | B73 reference; tissue development |
| Arabidopsis thaliana | AtGenExpress / eFP | bar.utoronto.ca | 100+ | Gold standard; hundreds of conditions |
| Arabidopsis | TRAVA | travadb.org | 30+ | Tissue-specific expression variance |
| Wheat (Triticum aestivum) | expVIP | wheat-expression.com | 20+ | Hexaploid; developmental stage series |
| Soybean (Glycine max) | SoyBase / eFP | bar.utoronto.ca/efpsoybean | 20+ | Nodule and root emphasis |
| Tomato (Solanum lycopersicum) | TomExpress / eFP | bar.utoronto.ca/efptomato | 20+ | Fruit development stages |
| Barley (Hordeum vulgare) | BARLEX | barlex.barleysequence.org | 15+ | Tissue panel |
| Cotton (Gossypium) | CottonGen | cottongen.org | 15+ | Fiber development emphasis |
| Medicago truncatula | Medicago eFP | bar.utoronto.ca/efpmedicago | 15+ | Nodule time course |
| Brassica rapa | B rapa eFP | bar.utoronto.ca/efpbrassica | 15+ | Morphotype comparison |

### Using Expression Atlases for eQTL

1. **Validate tissue specificity**: If an eGene is called "leaf-specific eQTL" but the
   gene is actually only expressed in roots → question the call.
2. **Identify top candidate genes**: In an eQTL peak region, prioritize genes with:
   - High expression in the sampled tissue
   - Expression variation across genotypes (not constitutively expressed)
   - Known function related to the trait of interest
3. **Cross-species inference**: If your species lacks an atlas, use Arabidopsis and rice
   orthologs. Conserved expression patterns across species increase confidence.

## Non-Model Species eQTL Strategy

When working with species not listed above:

1. **Cis-eQTL only**: Without a well-annotated genome and expression atlas, trans-eQTL
   interpretation is nearly impossible. Focus on cis.
2. **Homology-based annotation**: BLAST significant eGenes against Arabidopsis/rice
   proteomes. Prioritize eGenes with strong E-value hits (< 1e-10) to known plant genes.
3. **RNA-seq as pseudo-atlas**: If you have expression data from a single tissue only,
   acknowledge this limitation. A gene may be an eGene in your sampled tissue but not
   in others — this is expected, not a sign of false positives.
4. **Candidate gene approach**: If n < 50 and genome is fragmented (scaffold-level),
   restrict analysis to gene regions with known functional annotations.

## Statistical Power Guidelines

### cis-eQTL Power by Sample Size (Plant Studies)

| Sample Size | Detectable Effect (R²) | Expected eGenes (FDR < 0.05) | Notes |
|-------------|------------------------|------------------------------|-------|
| 30-50 | > 0.30 | 0-50 | Only large-effect cis-eQTL |
| 50-100 | > 0.15 | 100-1000 | Typical plant study range |
| 100-200 | > 0.08 | 500-5000 | Good cis coverage |
| 200-500 | > 0.05 | 2000-15000 | Trans-eQTL becomes feasible |
| 500+ | > 0.03 | 5000-30000 | Large crop panel power |

### trans-eQTL Power

Trans-eQTL detection requires approximately 5-10x larger sample sizes than cis-eQTL for
equivalent effect sizes. Sample size < 200 is generally insufficient for trans-eQTL
with acceptable FDR.

## eQTL Hotspot Interpretation

An eQTL "hotspot" is a genomic locus associated with many distant genes:

1. **True regulatory hotspots**: Often contain master transcription factors (e.g., MYB,
   WRKY, NAC family members in plants). These are biologically meaningful.
2. **Technical artifacts**: Can arise from batch effects, population stratification,
   or cross-homeolog mapping in polyploids.
3. **Validation**: Check whether hotspot genes are enriched for known TF binding motifs.
   Cross-reference with ChIP-seq or DAP-seq data if available for the species.

## Environmental and Developmental Context

Unlike human eQTL studies that can control environment, plant studies must account for:

- **Field vs controlled environment**: Field samples have much higher expression variance.
  eQTL effects tend to be diluted in field conditions but are more ecologically relevant.
- **Developmental stage**: A gene may only show eQTL effects at specific developmental
  stages. Time-series eQTL designs are becoming more common in plant research.
- **Diurnal/circadian effects**: Sampling time matters. Photosynthesis-related genes
  show strong diurnal expression patterns that can mask or mimic eQTL effects.
- **Stress-responsive eQTL**: Many plant eQTLs are only detectable under stress (drought,
  salt, pathogen). Consider this when interpreting "negative" results from single-condition
  studies.

## Key References

- GTEx Consortium approach adapted for crops: use PEER factors + MatrixEQTL + mashR
- Plant eQTL review: Druka et al., "Expression quantitative trait loci in plants"
- Polyploid eQTL: Zhang et al., "eQTL mapping in polyploid crops"
- Multi-tissue: The GTEx multi-tissue eQTL framework applied to plant atlases
