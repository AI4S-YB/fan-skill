# Time-Series Expression: Plant-Specific Considerations

## Why Plant Time-Series Differs from Mammalian

| Aspect | Mammalian / Generic | Plant Specific |
|--------|--------------------|---------------|
| Sampling | Controlled lab environment | Field diurnal variation, temperature fluctuation |
| Temporal scale | Minutes to hours (cell lines) | Days to months (developmental) |
| Circadian regulation | Limited to clock neurons | 30-50% of all genes oscillate |
| Replicates | Biological replicates only | Biological + plot/field technical replicates |
| Tissue complexity | Organ-specific cell types | Plastid + nuclear gene coordination |
| Polyploidy | Rare | Common; homeolog expression patterns |
| Stress response | Acute and resolved | Chronic stress adaptation common |

## Plant Developmental Stages

### Model species reference stages

| Species | Stage system | Key transitions | Data source |
|---------|-------------|-----------------|-------------|
| Arabidopsis | Boyes et al. 2001 | Vegetative → Bolting → Flowering → Silique | AtGenExpress, TraVA |
| Rice (japonica) | Zadoks / BBCH scale | Seedling → Tillering → Booting → Heading → Grain fill | RiceXPro, eFP |
| Rice (indica) | BBCH | Same as japonica, 1-2 weeks faster | Rice Expression Database |
| Maize | Iowa State leaf grading | V3-V18 vegetative, R1-R6 reproductive | MaizeGDB, qTeller |
| Soybean | Fehr & Caviness 1971 | VE-Vn vegetative, R1-R8 reproductive | SoyBase, SoyKB |
| Wheat | Zadoks scale (0-99) | Tillering → Stem elongation → Booting → Anthesis → Milk → Dough | expVIP, Wheat Expression Atlas |
| Tomato | BBCH | Seedling → Flowering → Fruit set → Ripening | TomExpress, SGN |
| Cotton | Oosterhuis 1990 | Emergence → Squaring → Flowering → Boll development | CottonGen |

### Stage metadata for time-series design

When designing a time-series experiment, define:
1. **Sampling interval**: hours (stress), days (development), weeks (seasonal)
2. **Replicate structure**: per time-point biological replicates (3 minimum)
3. **Time-of-day**: for field samples, always record sampling time to account for diurnal effects
4. **Developmental synchrony**: tag/select plants at the same developmental stage (not just the same chronological age)

## Tissue Atlases for Major Crops

### Public expression atlases

| Resource | Species | Tissues / Stages | URL |
|----------|---------|-----------------|-----|
| eFP Browser (BAR) | Arabidopsis, rice, maize, soybean, poplar, Medicago, potato, tomato, barley | 10-100+ per species | bar.utoronto.ca |
| RiceXPro | Rice (japonica) | Field transcriptome: 15 tissues, 7 developmental stages | ricexpro.dna.affrc.go.jp |
| MaizeGDB qTeller | Maize | 44 tissues, B73 reference | qteller.maizegdb.org |
| expVIP | Wheat, barley | 1,000+ RNA-seq samples | wheat-expression.com |
| TraVA | Arabidopsis | 79 tissues/stages, developmental series | travadb.org |
| TomExpress | Tomato | 100+ tissues, genotypes | tomexpress.versailles.inra.fr |
| PLAZA | 70+ plant species | Comparative expression | bioinformatics.psb.ugent.be/plaza |
| ATTED-II | Arabidopsis, rice, soybean, poplar | Co-expression networks | atted.jp |
| PlantGenIE | Poplar, spruce | Tissue atlas, developmental series | plantgenie.org |

### Using atlases for validation

1. Cross-reference cluster members against tissue atlas expression.
2. A cluster of "root-specific" genes in your leaf time course is a red flag.
3. Use atlas data to annotate clusters: "Cluster 3 resembles the maize leaf
   maturation gradient (qTeller leaf gradient data)."

## Circadian Gene Handling

### Detection

Circadian periodicity detection before clustering is essential for
plant time-series:

```r
# MetaCycle for circadian detection
library(MetaCycle)

write.table(data.frame(
  gene = rownames(expr_matrix),
  expr_matrix
), "metacycle_input.txt", sep = "\t", quote = FALSE, row.names = FALSE)

meta2d(
  infile = "metacycle_input.txt",
  filestyle = "txt",
  timepoints = time_points,
  outputFile = TRUE,
  outDir = ".",
  minper = 20,  # minimum period (hours)
  maxper = 28,  # maximum period (hours)
  ARSdefaultPer = 24,
  outIntegration = "both"
)
```

### Key decisions

| Scenario | Action |
|----------|--------|
| Dense time course, known circadian species | Run MetaCycle first. Separate circadian from non-circadian genes. Analyze separately. |
| Sparse time course (< 6 points) | Do NOT attempt periodicity detection. Flag in report: "insufficient temporal resolution for circadian analysis." |
| Field data with temperature fluctuation | Temperature-driven oscillation can mimic circadian. If temperature log available, run correlation of expression with temperature. |
| Constant-light or constant-condition | True circadian signal. Report period, phase, amplitude. |
| Light-dark cycle | Confounded: circadian + light response. Cannot separate without constant conditions. Acknowledge limitation. |

### Plant circadian genes (benchmark)

Known clock genes to validate your periodicity detection:
- Arabidopsis: CCA1, LHY, TOC1, PRR5/7/9, GI, ZTL, ELF3/4, LUX
- Rice: OsCCA1, OsLHY, OsTOC1, OsPRR1/37/73/95, OsGI, OsZTL
- Maize: ZmCCA1, ZmLHY, ZmTOC1, ZmPRR37/73, ZmGI

If your pipeline does NOT detect these as periodic (p < 0.05), either
the sampling scheme is inadequate or the algorithm parameters need tuning.

## Polyploidy and Homeolog Expression

### Allopolyploids (wheat, cotton, canola, tobacco)
- Subgenome-specific expression bias is common and biologically meaningful.
- Run clustering separately per subgenome, then overlay.
- Check if homeolog pairs fall into the same or different temporal clusters.
  Different clusters = subfunctionalization in temporal expression.

### Autopolyploids (potato, alfalfa, sugarcane)
- Homeologs are not distinguishable by sequence.
- Expression is the sum/dosage across all copies.
- Higher expression variance is expected; wider SE ribbons on line plots.

## Stress Response Time Courses

| Stress | Peak response window | Key markers (Arabidopsis) |
|--------|---------------------|--------------------------|
| Drought | 2-6 hr (ABA-dependent), 24-48 hr (adaptation) | RD29A, RD29B, DREB2A, LEA |
| Salt | 1-6 hr (osmotic), 24-72 hr (ionic) | SOS1, NHX1, HKT1 |
| Cold | 1-4 hr (CBF regulon), 24 hr+ (COR genes) | CBF1/2/3, COR15A, COR47 |
| Heat | 0.5-2 hr (HSP), 6-24 hr (thermotolerance) | HSP70, HSP101, HSFA2 |
| Pathogen (PTI) | 0.5-1 hr (ROS burst), 2-6 hr (defense) | FRK1, WRKY33, PR1 |
| Pathogen (ETI) | 2-6 hr (HR), 12-24 hr (SAR) | PBS1, RPM1, NPR1 |
| Wounding | 0.5-2 hr (JA), 6-24 hr (PI) | JAZ, MYC2, PIN2 |
| Nutrient deficiency | 24 hr to days | NRT, PHT, IRT1, bHLH TFs |

### Sampling window advice
- Stress experiments should include a "pre-treatment" (t0) sample from
  each biological replicate as its own control.
- The most informative window is usually the first 24 hours.
  If budget limits time points, bias toward early time points.

## Tissue-Specific Expression in Plants

### Organs and their marker genes

| Organ/Tissue | Marker genes (Arabidopsis) | Crop orthologs |
|-------------|---------------------------|-----------------|
| Leaf mesophyll | RBCS, LHCB, CAB | Same in all plants |
| Root cap | SMB, BRN1/2 | Check root-specific in crop atlas |
| Root hair | RHS, EXP7 | -- |
| Vascular bundle | VND7, SND1 | -- |
| Shoot apical meristem | STM, WUS, CLV3 | Maize: KN1 |
| Flower meristem | LFY, AP1, FUL | Rice: RFL, OsMADS14 |
| Pollen | LAT52, LTP12 | Highly conserved across angiosperms |
| Endosperm | FIS2, MEA | Cereals: specific storage proteins |
| Aleurone | LTP, GAMYB | Barley/wheat aleurone-specific |
| Nodule | ENOD40, Lb | Legume-specific |

Use these as positive controls for tissue-specificity analysis:
a gene with high tau in leaf samples should include RBCS in the top decile.

## References for Further Reading

- Mockler TC et al. (2007) "The DIURNAL project: ..." -- circadian expression database
- Stelpflug SC et al. (2016) "An Expanded Maize Gene Expression Atlas..." -- maize B73 atlas
- Sato Y et al. (2013) "RiceXPro Version 3.0..." -- rice field transcriptome
- Borrill P et al. (2016) "expVIP: a Customizable RNA-seq Data Analysis and Visualization Platform" -- wheat
- Klepikova AV et al. (2016) "A high resolution map of the Arabidopsis thaliana developmental transcriptome..." -- TraVA
