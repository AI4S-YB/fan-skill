# Plant Pan-genome Specific Considerations

## Plant Genome Size and Ploidy

### Small Genomes (< 500 Mb)
- **Examples**: rice (400 Mb), Brachypodium (270 Mb), Setaria (500 Mb)
- **Approach**: Full chromosome-level pan-genome feasible for 50+ accessions
- **Tool**: PGGB with `-p 95 -s 100k` for high-resolution graph

### Medium Genomes (500 Mb - 2 Gb)
- **Examples**: soybean (1 Gb), maize (2.3 Gb), sorghum (730 Mb)
- **Approach**: Chromosome-by-chromosome or split by synteny blocks
- **Tool**: Minigraph per chromosome, then merge

### Large/Complex Genomes (> 2 Gb)
- **Examples**: wheat (17 Gb, hexaploid), sugarcane (10 Gb, octoploid)
- **Approach**: Subgenome separation first (e.g., IWGSC RefSeq for wheat)
- **Tool**: Cactus with subgenome alignment, or limited to one subgenome

## Polyploidy in Pan-genomics

### Allopolyploids (wheat, cotton, canola, peanut)
- Separate subgenomes before pan-genome construction
- Each subgenome has its own core/variable gene classification
- Homeolog-specific PAV: gene present in A subgenome but absent in B subgenome
- PAV between subgenomes is expected (subgenome fractionation)

### Autopolyploids (potato, alfalfa, sugarcane)
- Haplotype-resolved assemblies needed for accurate PAV
- Allele dosage complicates "presence/absence" calling
- Use PanGenie with population-level k-mer frequencies

## PAV and Functional Impact

### Disease Resistance Genes (NBS-LRR)
- Largest PAV category in most plant species
- Rice: ~30% of NBS-LRR genes show PAV between indica and japonica
- Maize: >40% of NBS-LRR genes are variable across inbred lines
- **Implication**: PAV in R genes directly impacts disease resistance phenotyping

### Secondary Metabolism
- Terpene synthases (TPS), cytochrome P450s: highly variable
- Often arranged in tandem arrays that expand/contract between accessions
- PAV in these clusters explains metabolic diversity

### Transposable Elements (TEs)
- TE insertions near genes can cause "presence" of a gene variant
- Plant genomes: TEs are 50-85% of genome content
- PAV caused by TE insertion needs orthogonal validation (PCR)

## Computational Strategy

### Memory Optimization for Plant Genomes

For large plant genomes, use strategies to reduce memory:

1. **Chromosome splitting**: Process each chromosome independently
2. **K-mer pre-filtering**: Use 21-mer coverage to identify conserved regions
3. **Anchor-based alignment**: Use conserved gene anchors from BUSCO

### Recommended Tools by Use Case

| Use Case | Tool | Reason |
|----------|------|--------|
| 10-20 genomes, quick survey | Minigraph | Fast, good SV detection |
| 20-50 genomes, full variation | PGGB | Complete variation graph |
| 50+ genomes, complex ploidy | Cactus | Scalable progressive alignment |
| Resequencing-based PAV | PanGenie | K-mer based, no assembly needed |
| SV genotyping | vg giraffe | Read mapping to pan-genome graph |

## Reference Databases for Plant Pan-genomes

- **Rice**: RPAN (http://cgm.sjtu.edu.cn/3kricedb/), Rice Pan-genome Browser
- **Soybean**: SoyBase (https://soybase.org/), PanSoy
- **Maize**: MaizeGDB (https://maizegdb.org/), NAM pan-genome
- **Wheat**: Wheat Pan-genome (10+ Genomes project), http://wheat.cau.edu.cn/
- **Tomato**: Sol Genomics Network (https://solgenomics.net/)
- **Brassica**: BRAD (http://brassicadb.cn/)
- **Cotton**: CottonGen (https://cottongen.org/)
