# Plant GWAS Specific Considerations

## Why Plant GWAS Differs from Human GWAS

| Aspect | Human GWAS | Plant GWAS |
|--------|-----------|------------|
| Population structure | Usually subtle (λ≈1.0-1.05) | Often strong (λ>1.1 normal for inbreds) |
| Kinship | Background relatedness | Known pedigree + population structure |
| LD | Long-range (EUR: ~50kb) | Variable: rice ~200kb, maize <1kb |
| Multiple testing | 5e-8 standard | May relax based on marker count |
| Replication | Independent cohorts | Multi-environment trials |
| Polyploidy | Rare | Common (wheat 6x, potato 4x, cotton 4x) |
| Candidate validation | GWAS → experimental follow-up | GWAS → breeding population validation |

## Species-Specific GWAS Recommendations

### Self-pollinated crops (rice, soybean, wheat, barley, tomato)
- **Always include K matrix** — population structure from breeding history is inherent
- Q (PCA) + K (kinship) model is the standard
- Consider subpopulation analysis if PCA reveals clear clusters
- HWE filtering: skip or use relaxed threshold

### Cross-pollinated crops (maize, sunflower, rye, alfalfa)
- LD decays fast → need high-density markers for fine mapping
- PCA alone may be sufficient for population structure
- NAM/MAGIC populations dramatically increase power
- HWE filtering: appropriate (random mating assumption holds)

### Polyploid crops
- **Allopolyploids** (wheat, cotton, canola): analyze subgenomes separately
- **Autopolyploids** (potato, alfalfa, sugarcane): require dosage-aware methods
- Marker naming conventions often encode subgenome information
- Power is diluted across homeologous loci

## Multi-Environment Trials (MET)

The standard in plant breeding. Always check if phenotype data has an environment column.

### Per-environment GWAS approach
1. Run GWAS within each environment
2. Compare significant SNPs across environments
3. "Stable QTL" = significant in ≥2 environments
4. "Environment-specific QTL" = significant in only 1 environment

### BLUP approach
1. Fit mixed model: trait ~ (1|genotype) + (1|environment)
2. Extract BLUPs per genotype
3. Run single GWAS on BLUPs

Which to choose:
- BLUP: when environments represent random samples of target population of environments
- Per-environment: when environments are structured treatments (e.g., drought vs. control)

## Centromere and Pericentromeric Regions

Many plant genomes have large repetitive centromeric regions where:
- Recombination is suppressed → long-range LD
- GWAS signals may be spurious
- Candidate gene identification is harder

**Recommendation**: Cross-reference significant SNPs with centromere coordinates from the reference genome annotation. Flag any peak within ±5Mb of the centromere.

## Genomic Selection Crossover

GWAS results can inform genomic selection:
- Significant SNPs can be used as fixed effects in GS models
- GWAS-identified QTL regions can be weighted in prediction models
- See `bio-plant-genomic-selection` Skill for the GS workflow
