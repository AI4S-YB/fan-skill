# Plant CRISPR/Cas9 Specific Considerations

## Plant Genome Editing Challenges

### Polyploidy
Plant polyploidy creates unique CRISPR design challenges:
- **Allopolyploids** (wheat, cotton, canola): Highly similar homeologous gene copies
  - One sgRNA often targets all subgenome copies (high conservation in coding regions)
  - Subgenome-specific editing requires targeting divergent regions (UTRs, introns)
  - Advantages: Can knock out all homeologs simultaneously
  - Disadvantages: Cannot study subgenome-specific function without specific sgRNAs

### Repetitive Genomes
Plant genomes are rich in repetitive sequences (50-85%):
- sgRNAs in repetitive regions have many off-targets
- Use RepeatMasker to filter sgRNAs in repeat regions
- Transposable element-derived sequences should be avoided for sgRNA design

### GC Content Variation
Plant genomes show extreme GC content variation:
- Grasses (maize, wheat): relatively uniform ~45% GC
- Some dicots show strong GC gradients along chromosomes
- GC-rich regions may have fewer NGG PAMs -> use Cas9-NG or Cas12a

## Transformation Methods by Species

| Species | Method | Tissue | Time to T0 | Notes |
|---------|--------|--------|------------|-------|
| Rice (japonica) | Agrobacterium | Callus from mature embryo | 3-4 months | High efficiency |
| Rice (indica) | Agrobacterium | Callus from mature embryo | 4-6 months | Lower efficiency than japonica |
| Maize (B73) | Agrobacterium | Immature embryo | 4-6 months | Genotype-dependent |
| Maize (Hi-II) | Agrobacterium | Immature embryo | 3-4 months | High transformation efficiency |
| Wheat | Particle bombardment | Immature embryo | 4-6 months | Genotype-dependent |
| Barley | Agrobacterium | Immature embryo | 3-4 months | Golden Promise most transformable |
| Soybean | Agrobacterium | Cotyledonary node | 3-5 months | Genotype-dependent |
| Tomato | Agrobacterium | Cotyledon | 2-3 months | Micro-Tom or M82 preferred |
| Potato | Agrobacterium | Internode/leaf | 3-4 months | Tetraploid complicates editing |
| Cotton | Agrobacterium | Hypocotyl | 5-7 months | Long regeneration time |
| Arabidopsis | Floral dip | Flower | 2-3 months | No tissue culture needed |
| Brassica napus | Agrobacterium | Hypocotyl | 4-6 months | Polyploid (AACC genome) |

## Cas9 Delivery Systems for Plants

### DNA-Based (T-DNA)
- **Binary vector**: 35S:Cas9 + U6/U3:sgRNA in one T-DNA
- **Advantage**: Simple, stable expression
- **Disadvantage**: Continuous Cas9 expression -> more off-targets, transgene integration

### RNP (Ribonucleoprotein) Complexes
- Cas9 protein + in vitro transcribed sgRNA delivered via particle bombardment or PEG-mediated protoplast transfection
- **Advantage**: Transient, no transgene integration -> transgene-free editing
- **Disadvantage**: Lower efficiency, requires protoplast regeneration

### Virus-Based (VIGE - Virus-Induced Genome Editing)
- TRV, BSMV, or CPMV vectors delivering sgRNA into Cas9-transgenic plants
- **Advantage**: Bypasses tissue culture for sgRNA delivery
- **Disadvantage**: Requires Cas9-transgenic line, virus host range limits

## sgRNA Expression Systems in Plants

### RNA Polymerase III Promoters
| Promoter | Species | Expression Level |
|----------|---------|-----------------|
| OsU3 | Rice | High |
| OsU6a/OsU6b | Rice | High |
| TaU6 | Wheat | High |
| AtU6-26 | Arabidopsis | High |
| GmU6 | Soybean | Moderate |

### RNA Polymerase II Promoters (for multiplexing)
- **Csy4 system**: Csy4 ribonuclease cleaves polycistronic transcript
- **tRNA system**: Endogenous tRNA processing releases individual sgRNAs
- **Ribozyme system**: Self-cleaving hammerhead ribozymes

## Plant-Specific Editing Outcomes

### NHEJ Repair Spectrum in Plants
Plants show a distinct indel spectrum compared to animals:
- **1 bp insertion**: Most common (30-40% of edits) — typically +T or +A
- **1-3 bp deletion**: Common (25-35%)
- **Larger deletions**: Less common (10-20%)
- **Complex indels**: Deletion + insertion at same site

### HDR in Plants
HDR is very inefficient in plants (<1-5%):
- Strategies to improve:
  1. **Geminivirus replicon**: Amplify donor template in planta -> up to 10x improvement
  2. **NHEJ inhibition**: Scr7, Nu7441 chemical inhibitors, or ku70 mutants
  3. **Cell cycle synchronization**: HDR is active in S/G2 phase
  4. **Cas9 fusion**: Cas9-ctIP or Cas9-MRE11 fusions recruit repair factors
- Alternative to HDR: **Prime Editing** — no donor template needed, can make precise small edits (1-40 bp)

### Base Editing in Plants
- **CBE (Cytosine Base Editor)**: C->T conversion, editing window positions 4-8
  - Plant-optimized: rAPOBEC1-Cas9n-UGI
  - A3A-Cas9n-UGI: wider editing window (positions 2-13)
- **ABE (Adenine Base Editor)**: A->G conversion, editing window positions 4-7
  - Plant-optimized: TadA8e-Cas9n
  - PABE-7: evolved for plant use
- **Key for plants**: Need to introduce STOP codons (CAA/CAG/CGA->TAA/TAG/TGA)

## Recommended Databases and Resources

- **CRISPOR** (http://crispor.tefor.net/): sgRNA design with plant genomes
- **CRISPR-P** (http://crispr.hzau.edu.cn/CRISPR2/): Plant-focused CRISPR design
- **Cas-OFFinder** (http://www.rgenome.net/cas-offinder/): Off-target search
- **Plant genome databases**: Phytozome (https://phytozome-next.jgi.doe.gov/), Ensembl Plants
- **Addgene**: Plant CRISPR plasmids repository (https://www.addgene.org/)
- **CRISPR-GE** (http://skl.scau.edu.cn/): Toolkit for plant genome editing analysis
