# Plant Population Genetics Specific Considerations

## Inbreeding Coefficients

Self-pollinating crops have high inbreeding coefficients (F ≈ 1.0 for highly inbred lines). This affects:
- Expected heterozygosity: much lower than HWE expectations
- PCA: inbred individuals appear more "extreme" on PC axes
- Fst: inbreeding inflates Fst estimates (use Weir & Cockerham 1984 method)

## Polyploidy in Population Genetics

- Allopolyploids: run PCA per subgenome
- Autopolyploids: allele dosage coding matters
- Paleopolyploids (soybean, maize): treat as diploid for most analyses

## Relatedness in Breeding Populations

Modern breeding populations contain complex relatedness:
- Full-sibs, half-sibs in biparental crosses
- Shared parents across multiple crosses
- This creates structure that is real (not confounding) but must be modeled correctly
