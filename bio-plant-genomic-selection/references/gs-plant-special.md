# Plant Genomic Selection Specific Considerations

## Inbred vs Outcross Species

### Self-pollinated crops (rice, wheat, soybean, barley)
- Additive effects dominate → GBLUP/rrBLUP near optimal
- Non-additive variance is small in homozygous lines
- G×E is the main challenge, not model complexity
- Long-term GS: monitor loss of genetic variance

### Cross-pollinated crops (maize, sunflower, rye)
- Hybrid prediction: GCA vs SCA estimation
- Dominance and epistasis may matter for hybrid performance
- Training population design: include both parents and hybrids
- Heterotic groups must be modeled explicitly

## Multi-Environment GS

The standard in plant breeding. Options:

1. **Within-environment GS**: train and predict within each environment
2. **Across-environment GS**: train on all environments, predict new environment
3. **G×E GS**: explicit G×E covariance modeling (recommended for >3 environments)

For variety recommendation: across-environment GS with G×E
For genomic selection within a breeding program: within-environment often sufficient

## Training Population Optimization

- **Diversity**: training set should represent the breeding germplasm
- **Relatedness**: prediction accuracy drops when training and selection candidates are unrelated
- **Update frequency**: retrain every 2-3 breeding cycles as allele frequencies shift

## Genomic Selection vs Marker-Assisted Selection

- MAS: select on few large-effect QTL
- GS: select on genome-wide markers (all small effects)
- In practice: GS + major QTL as fixed effects = best approach for most traits

## Long-term Genetic Gain

Factors limiting long-term GS response:
- Loss of genetic variance (Bulmer effect)
- Recombination breaking favorable haplotypes
- G×E interaction in target environments

Mitigation strategies:
- Optimal contribution selection (balance gain vs diversity)
- Weighted GS (up-weight rare favorable alleles)
- Rapid cycling (speed breeding + GS)
