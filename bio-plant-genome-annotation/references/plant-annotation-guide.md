# Plant Genome Annotation Guide

## Why Plant Annotation Differs from Animal Annotation

| Aspect | Animal Genomes | Plant Genomes |
|--------|---------------|---------------|
| Repeat content | Moderate (30-50%) | High to very high (50-90%) |
| Gene density | ~10-12 genes/Kb | ~1-5 genes/Kb (varies by genome size) |
| Intron size | Generally small | Can be very large (especially grasses) |
| Polyploidy | Rare | Common (crop species frequently polyploid) |
| RNA editing | Rare in nucleus | Organellar RNA editing common |
| UTR length | Well characterized | Highly variable, poorly annotated |
| Pseudogenes | Common | Common, often from polyploidy |

## Species-Specific Annotation Recommendations

### Small-genome dicots (Arabidopsis thaliana, ~135Mb)
- Repeat masking is fast; use sensitive mode
- RNA-seq from multiple tissues is readily available
- Serves as the gold standard for plant annotation

### Medium-genome monocots (Oryza sativa, ~430Mb)
- Repeat masking takes moderate time
- Use Poales-specific BUSCO lineage
- High-quality reference annotations exist for comparison

### Large-genome monocots (Zea mays, ~2.4Gb; Hordeum vulgare, ~5Gb)
- >80% repeats; masking is the bottleneck
- Partition genome by chromosome for parallel processing
- EDTA pipeline may be preferred over RepeatModeler

### Gigantic genomes (Triticum aestivum, ~17Gb; Pinus taeda, ~22Gb)
- Must run masking and annotation per chromosome/scaffold
- Partition strategy essential for all steps
- Use subgenome-specific evidence for polyploids

### Polyploid crops
- **Allotetraploids** (Brassica napus, Gossypium hirsutum): Annotate subgenomes separately
- **Allohexaploids** (Triticum aestivum): 3 subgenomes; gene count ~3x diploid relatives
- **Autopolyploids** (Solanum tuberosum, Medicago sativa): Haplotype-resolved annotation needed

## Repeat Masking Strategy

### De novo vs Library-based
- **De novo** (RepeatModeler): Essential for non-model species; finds species-specific TE families
- **Library-based** (RepeatMasker with known repeats): Faster; good when close relative already characterized

### Plant-Specific Tools
- **EDTA**: Integrated pipeline specifically designed for plant genomes; handles LTR-RTs better than RepeatModeler alone
- **LTR_FINDER_parallel**: Scalable LTR retrotransposon discovery
- **LTR_retriever**: Post-processing for LTR candidates

## Gene Prediction Evidence Hierarchy

1. **Full-length cDNA/RNA-seq** (best): Precise exon-intron boundaries
2. **Iso-Seq/PacBio long reads**: Full-length isoforms without assembly
3. **Short-read RNA-seq** (multiple tissues): Good coverage, assembly needed
4. **Protein homology** (OrthoDB, UniProt/Swiss-Prot): Captures conserved genes only
5. **Ab initio** (no evidence): Least reliable; last resort

## Functional Annotation Best Practices

1. Always use both EggNOG-mapper and InterProScan — they are complementary
2. Use plant-specific taxonomic scopes for orthology assignment
3. Run KEGG pathway mapping with plant-specific KEGG databases (kegg_plant)
4. For crop species, PlantCyc provides better metabolic pathway coverage than KEGG alone
5. GO enrichment analysis requires plant-specific GO annotation database

## QC and Validation

### Internal Consistency Checks
- Gene count within 20% of closest annotated relative
- Average CDS length ~1000-1500bp
- Average exon count 4-6 for dicots, 4-8 for monocots
- < 5% of genes as single-exon (unless validated by RNA-seq)

### Cross-validation
- Compare BUSCO scores with related species' published annotations
- Check synteny with close relatives (if chromosome-level assembly available)
- Validate gene models with independent RNA-seq datasets
