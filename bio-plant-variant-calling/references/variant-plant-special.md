# Plant Variant Calling — Unique Considerations

## Why Plant Variant Calling Differs from Human

| Aspect | Human | Plant |
|--------|-------|-------|
| Genome size | 3.2 Gb (compact) | 0.1-22 Gb (rice 400Mb to pine 22Gb) |
| Ploidy | Diploid | Diploid to dodecaploid (sugarcane 12x) |
| Repeat content | ~50% | 60-95% (wheat ~85%, maize ~85%) |
| Heterozygosity | ~0.1% | 0-30% (inbred <0.1%, outcross up to 30%) |
| Genomic resources | Excellent | Variable (well-annotated to draft scaffolds) |
| Population panels | 1000 Genomes, gnomAD | Species-specific, often unpublished |
| Organelle DNA | <1% of reads | 5-20% of total reads |
| LD | 5-50 kb | 1-500 kb (highly species-dependent) |

## Polyploid Calling Strategies

### The Polyploid Problem

In a polyploid genome, a single read from a given locus could originate from any of the N homeologous copies. This creates three problems:
1. **Mapping ambiguity**: reads may map to the wrong homeolog
2. **Dosage uncertainty**: at a heterozygous site, the read count does not directly translate to genotype
3. **False heterozygosity**: reads from two different homeologs align to the same position, creating apparent heterozygosity where none exists

### Allopolyploid Strategy (e.g., Wheat AABBDD, Cotton AADD)

1. **If subgenome references exist**: Align reads separately to A, B, and D genome references. This gives the most accurate mapping.
2. **If only whole-genome reference exists**: Align to the whole reference, then:
   - Filter reads by MAPQ > 30 to keep only confidently assigned reads
   - Use `-ploidy 2` per subgenomic region if you can partition by chromosome
   - Expect higher missing call rates

### Autopolyploid Strategy (e.g., Potato 4x, Alfalfa 4x)

1. **GATK with `-ploidy N`**: Simple but accuracy degrades above 4x
2. **freebayes with `--ploidy N`**: Often more accurate for autopolyploids; uses a more natural genotype likelihood model
3. **Specialized tools**:
   - `updog` (R package): Normalized UMI-based genotyping for polyploids
   - `polyRAD` (R package): Genotype calling from RAD-seq in polyploids
   - `fitPoly` (R package): Genotype calling from SNP array data in polyploids

### Quality Threshold Adjustments for Polyploids

| Filter | Diploid Threshold | Polyploid Adjustment |
|--------|-------------------|---------------------|
| Mapping quality (MQ) | > 40 | > 30 (expect lower due to multi-mapping) |
| Fisher strand bias (FS) | < 60 | < 80 (4x), < 100 (6x) |
| Quality by depth (QD) | > 2.0 | > 1.5 |
| Mapping rate expectation | > 90% | > 70% (4x), > 60% (6x) |
| Missing rate per site | < 20% | < 30% |

## Organelle Genome Variants

Plant resequencing captures abundant chloroplast and mitochondrial reads. Key considerations:

### Identification
Organelle contigs are typically named:
- Chloroplast: `ChrC`, `Pt`, `chloroplast`, `cpDNA`, `NC_*`
- Mitochondria: `ChrM`, `Mt`, `mitochondria`, `mtDNA`

### Variant Calling for Organelles
- **Separate from nuclear calling**: Use a separate BAM file containing only organelle-aligned reads
- **Depth expectations**: Organelle coverage can be 100-1000x. Set `--max-depth` high enough (e.g., 10000) to avoid excluding true variants
- **Heteroplasmy**: Organelles can have within-sample variant mixtures. Use `bcftools mpileup --per-sample-mF` to detect heteroplasmic sites
- **Ploidy**: Treat organelles as haploid (`-ploidy 1`). Organelles are generally uniparentally inherited in plants (maternal in angiosperms)

### Warning: Organelle Reads as Contaminants
If you are NOT interested in organelle variants, filter out organelle-mapped reads BEFORE variant calling. Otherwise:
- Ti/Tv ratio will be skewed
- Coverage statistics will be inflated
- Heterozygosity estimates will be biased

## Plant Ti/Tv Expectations

Unlike humans where Ti/Tv ~2.0-2.5 is a well-established benchmark, plants show more variation:

| Plant Group | Expected Genome-Wide Ti/Tv | Exonic Ti/Tv |
|-------------|---------------------------|--------------|
| Grasses (rice, maize, sorghum) | 1.8-2.2 | 2.2-2.8 |
| Legumes (soybean, common bean) | 1.7-2.2 | 2.0-2.6 |
| Brassicaceae (Arabidopsis, Brassica) | 1.9-2.4 | 2.3-2.9 |
| Solanaceae (tomato, potato) | 1.8-2.3 | 2.1-2.7 |
| Conifers (pine, spruce) | 1.5-2.0 | 1.8-2.3 |
| Polyploids (wheat, cotton) | 1.6-2.0 | 2.0-2.5 |

Reasons for lower Ti/Tv in plants:
1. Higher repeat content causes more false positive transversions
2. Larger effective population sizes in outcrossing species maintain more transversions
3. GC-biased gene conversion in grasses elevates transitions specifically
4. Polyploidy introduces additional apparent transversion signals from homeologous variation

### Interpreting Ti/Tv
- **< 1.5**: Strong indicator of excessive false positives. Investigate alignment, raw data quality, or filter parameters.
- **1.5-1.8**: Acceptable for polyploid or repeat-rich genomes. Check that exonic Ti/Tv is higher.
- **1.8-2.5**: Normal range for most plant species.
- **> 2.5**: Rare; may indicate overly aggressive false positive filtering. Check that real biological variants are not being discarded.

## Mapping-by-Sequencing (Mutant Identification)

A common plant-specific use case not covered by standard variant callers:

### Strategy
1. Cross mutant (in a given background) to a diverged accession
2. Pool F2 individuals showing the mutant phenotype
3. Sequence the pool
4. Call variants between the pool and the reference
5. Look for regions where only the mutant parent allele is present (homozygous mutant allele, no wild-type allele)

### Recommended Tools
- **SHOREmap**: Classical mapping-by-sequencing for Arabidopsis
- **MutMap**: For rice and other crops with a reference (EMS mutagenesis)
- **QTL-seq**: Bulked segregant analysis with extreme phenotype pools
- **NGSEP**: For non-model species with draft genomes

### Allele Frequency Expectations
- In the mutant pool at the causal locus: ~100% mutant allele (0% wild-type)
- At unlinked regions: ~50% each allele
- SNP index (mutant allele frequency) plot across chromosomes: sharp peak at the causal locus

## Structural Variants in Plants

Plants have abundant structural variants (SVs) — much more than humans:
- Transposable element (TE) insertions/deletions are the dominant SV type
- Presence/absence variation (PAV): entire genes can be present in one accession and absent in another
- Copy number variation (CNV): resistance gene clusters commonly vary in copy number

This Skill focuses on SNPs/InDels. For SV detection:
- Short reads: Manta, DELLY, GRIDSS
- Long reads: Sniffles2, SVIM, cuteSV
- Population-level: Paragraph, GraphTyper
