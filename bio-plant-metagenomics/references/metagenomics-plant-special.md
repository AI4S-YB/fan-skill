# Plant Metagenomics Specific Considerations

## Plant Microbiome Compartments

### Rhizosphere (Root-Soil Interface)
- **Microbial density**: 10^8-10^10 cells/g soil
- **Community complexity**: High (1,000-10,000 species)
- **Challenges**: Plant DNA contamination 10-30%, soil humic acid inhibitors
- **Recommended sequencing depth**: >10 Gbp per sample
- **Key taxa**: Proteobacteria (Pseudomonas, Rhizobium), Actinobacteria (Streptomyces), Firmicutes (Bacillus)

### Endosphere (Root Interior)
- **Microbial density**: 10^4-10^7 cells/g tissue
- **Community complexity**: Low (10-100 species)
- **Challenges**: Plant DNA >80% of reads, low microbial biomass
- **Recommended sequencing depth**: >20 Gbp per sample
- **Key taxa**: Proteobacteria (Burkholderia, Herbaspirillum), Actinobacteria (Microbacterium)

### Phyllosphere (Leaf Surface)
- **Microbial density**: 10^3-10^6 cells/g tissue
- **Community complexity**: Very low (5-50 species)
- **Challenges**: Extremely low biomass, chloroplast DNA contamination, UV exposure
- **Recommended sequencing depth**: >15 Gbp per sample
- **Key taxa**: Proteobacteria (Methylobacterium, Sphingomonas), Bacteroidetes

## Host Read Removal Strategies

### Strategy A: Reference-Based Alignment (Preferred)

```bash
# Build plant genome index
bowtie2-build host_genome.fa host_index

# Align and filter
bowtie2 -x host_index -1 sample_R1.fq -2 sample_R2.fq \
  --un-conc-gz clean_sample.fq.gz \
  -p 16 > /dev/null
```

### Strategy B: Classification-Based Filtering

```bash
# If plant reference is incomplete
kraken2 --db plant_db --paired sample_R1.fq sample_R2.fq \
  --output kraken.out --report kraken.report

# Extract non-plant reads
extract_kraken_reads.py -k kraken.out \
  -s sample_R1.fq -s2 sample_R2.fq \
  -t 33090 --exclude \
  -o clean_sample_R1.fq -o2 clean_sample_R2.fq
```

### Efficiency Expectations

| Compartment | Before Removal | After Removal | Strategy |
|-------------|---------------|---------------|----------|
| Rhizosphere | 70-85% microbial | >95% microbial | Bowtie2 to host genome |
| Endosphere | 5-20% microbial | >80% microbial | Bowtie2 + Kraken2 double filter |
| Phyllosphere | 5-30% microbial | >75% microbial | Bowtie2 to chloroplast genome |
| Bulk Soil | >95% microbial | >98% microbial | Minimal removal needed |

## Plant-Specific Functional Genes

### Nitrogen Cycle
- **nifH**: Nitrogenase reductase — key for associative N-fixation in rhizosphere
- **amoA**: Ammonia monooxygenase — nitrification
- **nirK/nirS**: Nitrite reductase — denitrification
- **nosZ**: Nitrous oxide reductase — complete denitrification

### Plant Growth Promotion
- **acdS**: ACC deaminase — reduces plant ethylene stress
- **ipdC/iaaM**: IAA biosynthesis pathways
- **pqqC**: PQQ biosynthesis — phosphate solubilization
- **sidA**: Siderophore biosynthesis — iron acquisition

### Plant Pathogen/Defense
- **hrp/hrc**: Type III secretion system — pathogenic Pseudomonas/Xanthomonas
- **virA/virG**: Virulence regulon — Agrobacterium
- **chvE**: Chromosomal virulence — Rhizobiaceae

## Database Resources for Plant Microbiome

- **GREENGENES2**: Updated 16S + genome taxonomy (https://ftp.microbio.me/greengenes_release/)
- **SILVA**: Comprehensive rRNA database (https://www.arb-silva.de/)
- **GTDB**: Genome Taxonomy Database (https://gtdb.ecogenomic.org/)
- **PLANT_MICROBIOME_DB**: Collection of plant-associated metagenomes
- **MG-RAST**: Metagenomics Rapid Annotation (https://www.mg-rast.org/)
- **IMG/M**: Integrated Microbial Genomes & Microbiomes (https://img.jgi.doe.gov/)

## Recommended Sequencing Strategy

| Goal | Platform | Depth | Notes |
|------|----------|-------|-------|
| MAG recovery (soil) | Illumina NovaSeq 150PE | 20-30 Gbp | Deep sequencing needed for complex communities |
| MAG recovery (endosphere) | Illumina NovaSeq 150PE | 30-50 Gbp | Compensate for plant DNA |
| Long-read MAGs | ONT PromethION | 20 Gbp | Higher MAG contiguity |
| Functional survey | Illumina NovaSeq 100PE | 5-10 Gbp | Shorter reads adequate for gene catalog |
