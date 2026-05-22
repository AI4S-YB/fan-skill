# Plant RNA-seq Specific Considerations

## Polyploidy in RNA-seq

### Allopolyploids (wheat, cotton, canola, tobacco)
- **Homeolog expression bias**: A, B, D subgenomes in wheat often show unequal expression
- **Read mapping**: Map to combined reference with subgenome-specific references
- **Count assignment**: Use featureCounts with `-p` and careful GFF handling — homeologs need separate gene IDs
- **Bias analysis**: Compare DEG counts per subgenome — if one subgenome dominates DEGs, it may indicate functional dominance

### Autopolyploids (potato, alfalfa, sugarcane)
- **Allele dosage**: Expression differences may come from allele dosage, not true differential regulation
- **Reference choice**: Use haplotype-resolved reference if available

## Tissue/Organ Specificity

Plant gene expression is highly tissue-specific. Common contrasts:
- **Root vs Shoot**: Fundamental differentiation
- **Leaf developmental gradient**: Young → mature → senescent
- **Seed development**: Embryogenesis → maturation → desiccation
- **Stress response**: Drought, salt, heat, cold, pathogen
- **Circadian**: Time-of-day effects can be large — note sampling time

## Batch Effects

Common in plant RNA-seq:
- **Greenhouse vs Field**: Huge batch effect — always include as covariate
- **Sampling time of day**: Circadian genes vary >10-fold within 24h
- **Library prep date**: Standard batch effect, check PCA colored by prep date

## Low-Expression Gene Handling

Plants often have more tissue-specific genes than animals:
- ~30-50% of plant genes show tissue-specific expression
- Low-count filtering should use per-condition minimums (genes expressed in ANY condition should be retained)
- DESeq2's `results()` with `lfcThreshold=1` helps focus on biologically meaningful changes

## Reference Genome Quality

Critical for RNA-seq:
- **Good**: rice, maize, Arabidopsis — well-annotated, low missassembly
- **Adequate**: soybean, tomato, cotton — some gaps but usable
- **Challenging**: wheat (17Gb hexaploid), sugarcane (100+ chromosomes)
- For poorly annotated genomes: consider de novo transcriptome assembly (Trinity) as complement
