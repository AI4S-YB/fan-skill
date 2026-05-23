# Plant QTL Mapping — Unique Considerations

## Why Plant QTL Mapping Differs from Other Systems

| Aspect | Animal/Human QTL | Plant QTL |
|--------|-----------------|-----------|
| Population types | F2 intercross, backcross (rare) | F2, RIL, DH, BC, NAM, MAGIC — all common |
| Population permanence | Usually temporary | RIL/DH are permanent (immortal) |
| Replication | Limited | Multi-environment, multi-year common |
| Selfing | Rare | Common in self-pollinated crops |
| Marker density | High (SNP chips, WGS) | Variable (SSR → GBS → WGS) |
| Polyploidy | Rare | Common (wheat 6x, potato 4x, cotton 4x) |
| QTLxE analysis | Rare | Standard in plant breeding |
| Map construction | Physical map from reference | Genetic map from markers needed |

## Population-Specific Considerations

### F2 Populations
- Advantage: Quick to establish (one generation)
- Disadvantage: Cannot replicate, each individual is unique
- Genotype coding: AA (parent1), AB (heterozygote), BB (parent2)
- Can estimate additive + dominance effects
- Best for: Proof-of-concept, preliminary QTL mapping

### RIL Populations (Recombinant Inbred Lines)
- Advantage: Permanent, replicable, high recombination
- Disadvantage: Time-consuming (6-8 generations of selfing)
- Genotype coding: AA (parent1), BB (parent2) — no heterozygotes
- Only additive effects estimable (no dominance)
- Best for: High-resolution QTL mapping, multi-environment trials
- **Selfing vs Sib-mating**: Use `riself` for single-seed descent, `risib` for sib-mated RILs

### DH Populations (Doubled Haploids)
- Similar to RIL: homozygous, permanent
- Faster to produce (1-2 years vs 3-5 for RIL)
- No selection during production (unlike RIL where selection can occur)
- Available in: maize (most common), rice, wheat, barley, canola
- Not available in: soybean (low haploid induction), many legumes

### BC Populations (Backcross)
- BC1: F1 backcrossed to one parent
- Genotypes: AA (backcross parent homozygous) and AB (heterozygous)
- Only one parent's alleles are segregating
- Best for: Introgressing a specific trait from donor parent
- Less common in modern QTL studies

### MAGIC/NAM Populations (Multi-Parent)
- Multiple founders (4, 8, or more)
- Higher allelic diversity → detect more QTL alleles
- Higher mapping resolution (more historical recombination)
- Complex analysis: ASMap + R/qtl multi-parent functions
- Available in: rice (MAGIC), maize (NAM), wheat (MAGIC), Arabidopsis (MAGIC)

## Marker Type Considerations

### SSR/AFLP/RFLP (Low-Density, Legacy)
- 50-500 markers
- Low resolution (QTL CI: 15-30 cM)
- Use IM or CIM with few cofactors
- Map construction: R/qtl est.map (standard)

### SNP Chip (Medium-Density)
- 1,000-50,000 markers
- Medium resolution (QTL CI: 5-15 cM)
- Use CIM or MQM
- Map construction: ASMap or est.map

### GBS/WGS (High-Density)
- >50,000 markers
- High resolution (QTL CI: 1-5 cM)
- Use LepMap3 for map construction
- MQM or ICIM for QTL analysis
- High LD between adjacent markers → bin markers before analysis

## Multi-Environment QTL Analysis

Multi-environment data is a hallmark of plant QTL studies. The same RIL/DH population can be phenotyped in 3-10 environments (years x locations).

### Per-Environment QTL Comparison
1. Run QTL scan independently in each environment
2. Identify QTL in each environment
3. Classify:
   - **Stable QTL**: Detected in >= 2 environments, similar position and effect direction
   - **Environment-specific QTL**: Detected only in one environment
   - **Opposite-effect QTL**: Detected in >= 2 environments but opposite sign (rare, indicates strong QTLxE)

### Joint Multi-Environment Analysis
- Higher statistical power
- Estimates QTL main effect + QTLxE variance
- QTL with small QTLxE: good targets for marker-assisted selection
- QTL with large QTLxE: use for environment-specific breeding

### When to Use Which Strategy
- **Joint analysis**: When environments are random samples of target environment range
- **Per-environment**: When environments are structured treatments (drought vs irrigated, low N vs high N)
- **Both**: Run joint for power, per-environment for biological interpretation

## Segregation Distortion in Plants

Segregation distortion is common in plant populations, especially in wide crosses (interspecific or inter-subspecific).

### Causes
- Gametophytic selection (pollen competition)
- Zygotic selection (seedling lethality)
- Cytoplasmic-nuclear incompatibility
- Linkage to self-incompatibility loci

### Handling Strategies
1. **<10% distorted markers**: Ignore, keep all markers
2. **10-30% distorted markers**: Flag, run analysis with and without distorted markers
3. **>30% distorted markers**: Remove worst distorted markers (P < 0.001), re-run map
4. **Check for clustered distortion**: If distorted markers cluster in specific genomic regions → possible selection at linked viability loci

### Reporting
When publishing, report:
- Number and percentage of distorted markers
- Chromosomal distribution of distorted markers
- Whether removing distorted markers changes QTL results
- If distortion clusters overlap with QTL → discuss potential selection bias

## Centromeric and Pericentromeric QTL

In plants, QTL located near centromeres require special attention:
- Recombination is suppressed at centromeres → inflated LOD scores
- Confidence intervals are wider due to fewer recombination events
- Markers cluster at the genetic level but span large physical regions
- QTL peaks at centromeres may represent single large-effect genes OR artifacts

### Best Practices
1. Compare genetic vs physical map: is the peak at a recombination coldspot?
2. Check if the peak spans the centromere (broad flat LOD curve in that region)
3. If physical map available: calculate physical confidence interval separately from genetic
4. Report: "QTL on chrX spans centromeric region; fine-mapping requires additional recombinants"

## QTL Naming Conventions

### Standard Format
`q[TraitAbbreviation]-[Chromosome]-[Number]`

Examples:
- `qPH-1-1`: first plant height QTL on chromosome 1
- `qYLD-3-2`: second grain yield QTL on chromosome 3
- `qDTF-7-1`: first days-to-flowering QTL on chromosome 7

### Multi-Environment Extension
`q[Trait]-[Chromosome]-[Number].[Environment]`

Example:
- `qPH-1-1.E1`: plant height QTL specific to environment 1

## Integration with Other Analyses

### QTL + GWAS
- Co-localization of QTL and GWAS peaks strengthens candidate gene evidence
- QTL gives linkage-based support (coarse but high power)
- GWAS gives LD-based resolution (fine but lower power for rare alleles)

### QTL + RNA-seq
- eQTL mapping: treat gene expression as phenotype
- Candidate gene prioritization: genes within QTL interval that are differentially expressed

### QTL + Genomic Selection
- QTL results can be used as fixed effects in GS models
- Known large-effect QTL should be accounted for before GS
- See `bio-plant-genomic-selection` Skill
