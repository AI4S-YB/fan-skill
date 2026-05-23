# Plant Multi-Omics Specific Considerations

## Why Plant Multi-Omics Differs from Human/Animal Multi-Omics

| Aspect | Human/Animal | Plant |
|--------|-------------|-------|
| Sample types | Blood, tissue biopsies | Leaf, root, seed, fruit, fiber — fundamentally different metabolomes |
| Biological replicates | Clone/cell-line possible | Genetic heterogeneity even within inbred lines |
| Environmental control | Lab animals, controlled environments | Field conditions, seasonal variation, soil heterogeneity |
| Metabolome complexity | ~5,000 endogenous metabolites (HMDB) | >200,000 estimated (many unknown) specialized metabolites |
| Polyploidy | Rare | Common (wheat 6x, potato 4x, cotton 4x, canola 4x) |
| Genome annotation quality | High (human, mouse) | Variable — excellent (Arabidopsis, rice) to poor (orphan crops) |
| Temporal dynamics | Hours-days (circadian, fasting) | Minutes (stress response) to months (development, season) |
| Subcellular compartmentalization | Standard organelles | Adds chloroplasts, vacuoles, cell wall — metabolites partitioned |
| Tissue specialization | Moderate | Extreme — leaf vs. root vs. seed metabolomes near-orthogonal |
| Multi-omics studies | Clinical cohorts (n=100-1000s) | Typically n=10-100, often underpowered |

## Species-Specific Multi-Omics Recommendations

### Model species (Arabidopsis thaliana)
- Excellent genome annotation and TF databases
- Most published multi-omics reference datasets available
- Can use organism-specific OrgDb for enrichment (org.At.tair.db)
- Caveat: findings may not translate well to crops with different physiology

### Major cereals (rice, maize, wheat)
- **Rice**: Near-complete annotation, many public multi-omics datasets. Single reference genome.
- **Maize**: High genetic diversity (pan-genome). Consider using pan-genome references for RNA-seq mapping.
- **Wheat**: Hexaploid (A/B/D). Run subgenome-separated analyses. Homeolog co-expression networks often reveal polyploid-specific regulatory patterns.

### Legumes (soybean, common bean, chickpea, pea)
- Rich specialized metabolome: isoflavonoids, saponins, alkaloids
- Symbiotic nitrogen fixation adds root nodule as a unique tissue type
- Rhizosphere multi-omics (plant + microbiome) increasingly important

### Solanaceae (tomato, potato, pepper, eggplant)
- Fruit development and ripening are major research themes
- Potato: tuber metabolome dominated by starch and glycoalkaloids
- Tomato: well-characterized carotenoid and volatile pathways

### Oil crops (canola, sunflower, peanut, oil palm)
- Lipidomics integration requires specialized normalization (lipid class-specific)
- Seed development time-series multi-omics captures oil accumulation dynamics
- Canola (Brassica napus, allotetraploid A/C genomes)

### Fiber crops (cotton, flax, hemp)
- Fiber development stages have distinct multi-omics profiles
- Cotton: allotetraploid (A/D). Fiber initiation vs elongation vs secondary cell wall
- Cell wall metabolomics important for fiber quality

## Plant-Specific Data Integration Challenges

### 1. The "Unknown Metabolome" Problem

Plant metabolomes are dominated by unannotated features (often >50% of detected peaks). Strategies:

- **Keep unknown features**: Do not filter them out before MOFA2/DIABLO. They carry information.
- **Retrospective annotation**: After integration, prioritize unknown features that load heavily on biologically interesting factors for MS/MS-based identification.
- **In silico annotation**: Use tools like CANOPUS (SIRIUS) for compound class prediction.
- **Report unknowns honestly**: "Factor 3 is driven by 3 unknown features (m/z 285.076, rt 4.2 min; m/z 447.093, rt 5.1 min) and 2 known flavonoids."

### 2. Polyploid Subgenome Handling

- **Pre-processing**: Map RNA-seq reads to a combined reference that distinguishes subgenomes. Use tools like `HomeoRoq` or `polyCat` for homeolog-specific quantification.
- **Integration**: Run MOFA2 separately per subgenome OR label genes/features with subgenome prefix. Factor comparison across subgenomes reveals subgenome dominance patterns.
- **DIABLO with subgenome labels**: If subgenome origin can serve as a class label, DIABLO can identify features that distinguish subgenome contributions.

### 3. Tissue and Developmental Stage as Covariates

Plant multi-omics datasets often span multiple tissues or developmental stages. These introduce dominant sources of variation:

- **Do NOT pool tissues blindly**: A MOFA2 factor may primarily separate leaf from root — this is valid biology but may obscure treatment effects.
- **Stratify by tissue**: Run integration within each tissue, then compare factor structures.
- **Include as covariates**: Use `model_opts$covariates` in MOFA2 to regress out known tissue effects.

### 4. Circadian and Diurnal Effects

Many plant metabolic and transcriptomic processes follow circadian rhythms:

- **Record sampling time**: Essential metadata for plant multi-omics studies.
- **Time-stamped sampling**: If possible, sample all plants at the same time of day.
- **Post-hoc check**: After MOFA2, check if factor scores correlate with sampling time.

### 5. Stress Responses Dominate Multi-Omics Signal

Plant stress responses (drought, salt, pathogen, herbivory) produce massive multi-omics changes that can overwhelm more subtle signals:

- **Unstressed baseline**: Include adequate unstressed controls.
- **Stress-specific factors**: A MOFA2 factor explaining > 30% variance in transcriptome + metabolome is often a stress response.
- **Cross-stress comparison**: If multiple stress treatments, check for shared "general stress response" factors vs. stress-specific factors.

### 6. Soil and Rhizosphere Metadata

For field-grown plants or soil-based experiments:

- **Soil metadata matters**: pH, organic matter, nitrogen content, microbiome composition.
- **Multi-omics + environment**: Consider environment (soil properties) as an additional "omics" or metadata layer in integration.
- **Replication across environments**: Results from a single field may not generalize.

## Multi-Omics Integration for Plant Breeding

### Connecting integration results to breeding decisions

- **Stable multi-omics biomarkers**: Features consistently selected by DIABLO across environments (multi-environment MOFA/DIABLO) are candidate biomarkers.
- **Factor-based selection indices**: MOFA2 factor scores can be used as composite traits in genomic prediction models.
- **Omics-informed GS**: Weight genomic relationship matrices by feature importance scores from multi-omics integration.
- **Cross-trait integration**: Multi-omics can connect different traits (e.g., metabolome linked to yield components).

### From integration to functional validation

- **Prioritized gene list**: Top-weighted genes in MOFA2 factors that also correlate with a trait of interest are candidates for functional validation.
- **Network-guided mutagenesis**: Genes with high centrality in cross-omics networks are likely to have pleiotropic effects.
- **CRISPR targets**: Multi-omics prioritized genes are ideal CRISPR targets for crop improvement.

## Reference Multi-Omics Datasets

Key public plant multi-omics datasets to benchmark against:

- **Arabidopsis**: AtMORE (Multi-Omics Resource for Exploring gene function), 1001 Epigenomes
- **Rice**: RiceENCODE, RIS (Rice Information System) multi-omics modules
- **Maize**: MaizeGDB integrated datasets, NAM population multi-omics
- **Soybean**: SoyKB (Soybean Knowledge Base)
- **Tomato**: TomExpress + metabolomics atlas, Tomato Functional Genomics Database
- **Wheat**: expVIP (Expression Visualization and Integration Platform), WheatOmics
