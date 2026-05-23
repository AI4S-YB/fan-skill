# Plant Genotype Imputation: Plant-Specific Considerations

## Why Plant Imputation Differs from Human Imputation

| Aspect | Human Imputation | Plant Imputation |
|--------|-----------------|------------------|
| Reference panels | Large, diverse (1000G, TOPMed, HRC) | Species-specific, often limited or absent |
| LD extent | EUR ~50kb, AFR ~10kb | Highly variable: maize <1kb, rice ~200kb, wheat ~5Mb |
| Population structure | Subtle within continental groups | Often strong (breeding populations, landraces) |
| Ploidy | Diploid | Diploid to hexaploid (wheat) and higher |
| Inbreeding | Negligible | Common: selfing species (rice, soybean, wheat) |
| Genome size | ~3 Gb | Variable: rice 430 Mb to wheat 17 Gb |
| Tool support | All tools tested on humans | Limited testing; Beagle most used in plants |
| Phasing | SHAPEIT/Beagle | SHAPEIT may fail for exotic ploidy; Beagle preferred |

## When Imputation Is Useful in Plant Breeding

### Scenario 1: Merging Genotyping Platforms
Different breeding programs use different SNP arrays (e.g., 10K, 50K, 200K). Imputation to a common reference panel allows joint analysis.

**Action**: Impute all arrays to the highest-density platform or to a sequence-based reference.

### Scenario 2: Increasing Marker Density for GWAS
Low-density genotyping (e.g., GBS at 10K SNPs) produces wide GWAS peaks. Imputation increases resolution.

**Action**: Impute to sequence-level density, then run GWAS. Filter imputed SNPs at R² >= 0.5 for cautious GWAS.

### Scenario 3: Genomic Selection Input Preparation
GS models require all individuals to have the same SNP set. Different batches can be merged via imputation.

**Action**: Impute all batches, then unify to common variants with R² >= 0.3. Higher acceptance threshold for GS compared to GWAS.

### Scenario 4: Historical Data Integration
Older datasets with legacy markers (RFLP, SSR, low-density SNP) can be imputed to modern marker sets.

**Action**: Use LD information from a modern reference; flag that historical imputation accuracy is likely lower.

## Species-Specific Recommendations

### Rice (Oryza sativa) -- Self-pollinated, diploid
- **Reference panel**: 3K Rice Genomes Project (3024 accessions)
- **LD extent**: ~100-200 kb (indica), ~200-400 kb (japonica)
- **Recommended tool**: Beagle 5 (no reference needed for within-population) or Minimac4 (if using 3K reference)
- **Effective population size (Ne)**: ~100-200 for elite breeding lines
- **Special notes**: Highly structured (indica vs japonica vs aus); stratify by subspecies before imputation

### Maize (Zea mays) -- Cross-pollinated, diploid
- **Reference panel**: Maize HapMap3 (1218 lines), NAM founder lines
- **LD extent**: <1-2 kb (tropical), 10-100 kb (temperate)
- **Recommended tool**: Minimac4 or Beagle with large sample size
- **Effective population size (Ne)**: ~500-2000 for diverse populations
- **Special notes**: Fast LD decay means very high density reference is essential for accurate imputation

### Wheat (Triticum aestivum) -- Self-pollinated, allohexaploid (AABBDD)
- **Reference panel**: 10+ Wheat Genomes, Watkins landrace collection
- **LD extent**: 1-5 Mb within elite lines (extreme due to polyploidy + selfing)
- **Recommended tool**: Beagle (polyploid-aware) or split by subgenome + Beagle per subgenome
- **Effective population size (Ne)**: ~50-150 (strong domestication bottleneck)
- **Special notes**: ALWAYS split by subgenome (A/B/D) before imputation. Homeologous SNPs between subgenomes cause severe mis-phasing if not separated. Extended LD makes imputation relatively easy within each subgenome.

### Soybean (Glycine max) -- Self-pollinated, diploid (paleopolyploid)
- **Reference panel**: Soybean 302 accessions (USDA), wild soybean accessions
- **LD extent**: ~100-500 kb (elite), ~50 kb (wild)
- **Recommended tool**: Beagle or IMPUTE2 with panel
- **Effective population size (Ne)**: ~50-100 (extreme domestication bottleneck)
- **Special notes**: Consider including wild soybean (G. soja) in the reference for rare allele imputation

### Potato (Solanum tuberosum) -- Autotetraploid (2n=4x=48)
- **Reference panel**: Limited; SolCAP diversity panel
- **LD extent**: 1-5 Mb (clonal propagation + autopolyploidy)
- **Recommended tool**: Beagle (experimental polyploid mode) or polyRAD/UPdog for dosage calling first
- **Effective population size (Ne)**: ~100-300
- **Special notes**: Standard tools assume diploid. Dosage (0/1/2/3/4 copies) matters. Beagle's polyploid support is experimental. Verify results carefully.

### Cotton (Gossypium hirsutum) -- Allotetraploid (AADD)
- **Reference panel**: Cotton germplasm collections
- **LD extent**: ~100-500 kb
- **Recommended tool**: Beagle, split by subgenome (At/Dt)
- **Special notes**: Split At and Dt subgenomes as for wheat. Low polymorphism in Dt subgenome (single origin) means imputation accuracy is inherently lower for Dt.

### Non-Model Species
For species without a reference panel:
- **Only option**: Beagle 5 internal imputation
- **Minimum requirements**: >= 100 samples, >= 1000 SNPs
- **Expected accuracy**: Lower; DR2 typically 0.3-0.6
- **Honest reporting**: State that accuracy estimates (DR2) are based on internal cross-validation and may overestimate real accuracy
- **Mitigation**: Genotype a subset at higher density (e.g., 10% of samples by sequencing) for validation

## Ploidy-Aware Imputation Strategies

### Diploid Species
- All standard tools (Beagle, IMPUTE2, Minimac4) work well
- Standard VCF with GT field is sufficient
- Decision matrix: sample count + reference panel drives tool choice

### Allopolyploids (Wheat AABBDD, Cotton AADD, Canola AACC)
**Recommended workflow**:
1. Separate SNPs by subgenome using genome-specific marker annotation
2. Treat each subgenome as an independent diploid dataset
3. Impute each subgenome separately with Beagle
4. Re-merge after imputation

**Why this works**: Homeologous chromosomes in different subgenomes share sequence similarity but are genetically independent. Treating them as a single genome causes false LD signals between subgenomes.

### Autopolyploids (Potato 4x, Alfalfa 4x, Sugarcane 8x+)
**Challenge**: Standard VCF (0/1/2) does not encode dosage for polyploids.
**Solutions**:
1. Convert dosage to pseudo-diploid (0/1) -- loss of information but works with all tools
2. Use Beagle's experimental polyploid mode -- limited testing but most principled
3. Use polyRAD or fitPoly for dosage calling, then Beagle
4. Consider whether imputation adds value vs. analyzing the available markers directly

## Inbreeding and Its Effect on Imputation

### Self-Pollinated Crops
- Extended LD blocks make imputation easier within blocks
- But fewer recombination breakpoints mean phasing resolution is lower
- Lower Ne means fewer haplotype templates in the population
- **Practical implication**: Imputation accuracy is HIGH within LD blocks but phasing confidence is LOW at block boundaries

### Cross-Pollinated Crops
- Short LD blocks require more markers to capture haplotype diversity
- More recombination = better phasing resolution
- Higher Ne = more haplotype diversity
- **Practical implication**: Need higher starting marker density for adequate imputation (>= 10K SNPs minimum)

## Population Structure Management

Plant breeding populations often have strong structure:
- **Breeding program structure**: Elite x Exotic, Winter x Spring, etc.
- **Geographic structure**: Latitudinal adaptation, regional germplasm pools
- **Temporal structure**: Historical varieties vs. modern elite

**Before imputing**:
1. Run PCA on your target samples
2. If PC1 explains > 15% variance or clear visual clusters exist:
   - Stratify into subpopulations
   - Impute within each subpopulation separately
   - Merge after imputation
3. For reference-panel-based methods: ensure the reference covers all subpopulations

**After imputing**:
- Check that DR2 is consistent across subpopulations
- If one subpopulation has much lower DR2, its haplotypes were underrepresented in the reference/imputation template
- Consider re-running that subpopulation with a more matched reference

## Chunking Strategy for Plant Genomes

Plant genomes have vastly different sizes:

| Species | Genome Size | Recommended Chunk Size |
|---------|------------|----------------------|
| Arabidopsis | 135 Mb | 5 Mb |
| Rice | 430 Mb | 10-20 Mb |
| Soybean | 1.1 Gb | 20-30 Mb |
| Maize | 2.3 Gb | 30-50 Mb |
| Barley | 5.1 Gb | 50-100 Mb |
| Wheat | 17 Gb | 100-200 Mb per subgenome |

For very large genomes (wheat, barley), always:
- Split by chromosome/subgenome
- Process each chromosome chunk independently
- Use cluster computing with job arrays

## Validation Strategies for Plants

Without truth datasets, validate imputation by:

1. **Mask-and-impute**: Randomly mask 5-10% of known genotypes, impute, compare
2. **Cross-platform validation**: If samples were genotyped on multiple arrays, use one as truth
3. **Leave-one-out**: For Beagle internal imputation, the DR2 is an estimate of accuracy
4. **LD-based sanity check**: Post-imputation LD patterns should match expected LD decay
5. **Concordance with sequence data**: If any samples have sequence data, compare imputed calls
