# Plant GRN-Specific Considerations

## Why Plant GRN Differs from Animal GRN

| Aspect | Animal GRN | Plant GRN |
|--------|-----------|------------|
| TF repertoire | ~1600 TFs (human) | ~2000-3300 TFs (rice/maize), larger families |
| Dominant TF families | C2H2-ZnF, Homeobox, bHLH | MYB, WRKY, NAC, bHLH, AP2/ERF (plant-specific expansions) |
| Network architecture | Hierarchical, modular | Highly connected hubs, condition-specific rewiring |
| Motif databases | Extensive (JASPAR, CIS-BP) | Limited; PlantTFDB + Arabidopsis DAP-seq |
| Epigenetic regulation | CpG methylation | CG/CHG/CHH methylation, Polycomb, chromatin remodeling |
| Polyploidy | Rare | Common (wheat 6x, potato 4x, cotton 4x, canola 4x) |
| Experimental validation | ChIP-seq, Perturb-seq | DAP-seq, ChIP-seq, yeast one-hybrid (limited throughput) |
| Stress response | Less dominant in GRN | Central role — drought, pathogen, temperature stress rewiring |

## Major Plant TF Families

### MYB Superfamily

- The largest plant TF family (Arabidopsis: ~200 genes)
- Characterized by MYB DNA-binding domain (R1, R2, R3 repeats)
- R2R3-MYB: largest subfamily, regulates secondary metabolism, stress, development
- Regulation: ABA signaling, JA signaling, phenylpropanoid pathway
- GRN role: MYBs are consistently top hub TFs in most plant expression datasets

### WRKY Family

- Defined by WRKYGQK heptapeptide + zinc-finger motif
- Bind to W-box (TTGACC/T) in target gene promoters
- Central to biotic stress (pathogen defense, SA/JA crosstalk)
- Arabidopsis: 72 members; rice: >100 members
- GRN role: Major hubs in pathogen/stress response networks

### NAC Family (NAM/ATAF/CUC)

- Plant-specific — no animal orthologs
- Named after founding members: NAM, ATAF1/2, CUC2
- Bind to NACRS (NAC recognition sequence)
- Key roles: senescence, stress response, cell division, hormone signaling
- Rice: ~151 members; often top hubs in drought/stress GRNs

### bHLH (Basic Helix-Loop-Helix)

- Conserved across eukaryotes but massively expanded in plants
- Arabidopsis: 162 members; rice: 173
- Roles: JA signaling (MYC2), light signaling (PIFs), iron homeostasis
- GRN role: Medium hubs; often co-regulate with MYB and WRKY

### AP2/ERF Family

- Plant-specific DNA-binding domain (AP2 domain)
- Subfamilies: AP2, ERF, DREB, RAV
- DREB subfamily: binds DRE/CRT elements (drought/cold response)
- ERF subfamily: binds GCC-box (JA/ET pathogen response)
- GRN role: Stress-specific network hubs, particularly DREB1A/CBF3 in cold

### bZIP Family

- Basic leucine zipper domain, dimerization-dependent DNA binding
- Binds ABRE (ABA-responsive element)
- Key TFs: ABI5, AREB/ABF (ABA signaling core)
- GRN role: ABA response modules; stress hormone crosstalk hubs

### GRAS Family

- Named after founding members: GAI, RGA, SCR
- Includes DELLA proteins (GA signaling repressors)
- Roles: GA signaling, root development, stress adaptation
- GRN role: Hormone signaling crosstalk; developmental network hubs

### MADS-box Family

- Binds CArG-box motif [CC(A/T)6GG]
- Type I (SRF-like) and Type II (MEF2-like = MIKC-type)
- MIKC-type: floral organ identity (ABCDE model), flowering time control
- GRN role: Tissue-specific (floral) networks; often form multi-protein complexes

## Polyploid TF Network Considerations

### Allopolyploids (e.g., wheat AABBDD, cotton AADD)

- **Homeologous TF pairs**: Two or three copies of each TF gene from each subgenome
- **Expression partitioning**: Homeologs may have identical, diverged, or
  complementary expression patterns
- **Regulatory divergence**: HOX paralogs in wheat — A-genome copy regulates
  glume development, D-genome copy regulates flowering time
- **GRN analysis strategy**:
  1. Assign TFs to subgenomes using gene ID prefixes or synteny
  2. Build subgenome-specific GRNs AND a combined GRN
  3. Compare hub gene lists between subgenomes
  4. Identify "dominant homeologs" (hub in one subgenome but not the other)
  5. Check target overlap between homeologous TF pairs

### Autopolyploids (e.g., potato 4x, alfalfa 4x)

- No distinct subgenomes — allele dosage effects dominate
- GRN differences arise from allele-specific expression
- Recommend treating as diploid with increased expression variance

## Non-Model Species GRN Strategies

### TF Annotation

1. Primary: Sequence homology to Arabidopsis TFs (BLASTP/DIAMOND, pident >= 30%)
2. Secondary: PFAM domain scan for DNA-binding domains
3. Tertiary: iTAK (plant TF prediction pipeline) for automated plant TF annotation

### Network Validation

Without experimental data (ChIP-seq, DAP-seq), validate computationally:
1. Expression correlation enrichment (real edges should have higher correlation than random)
2. GO enrichment of target gene sets per TF
3. Cross-reference with Arabidopsis ortholog interactions (from STRING-DB, AtRegNet)
4. Motif enrichment in target promoters (if genome sequence available)

### Motif Analysis for Non-Model Species

- Extract promoter sequences (1-2 kb upstream of TSS) for each target gene
- Scan with Arabidopsis TF motifs (JASPAR plant collection)
- ENRICHEr-style enrichment: are motifs for TF-X enriched in TF-X's predicted targets?

## Network Quality Metrics for Plants

| Metric | Good | Borderline | Poor |
|--------|------|-----------|------|
| TFs identified | >1000 | 500-1000 | <500 |
| Edges (after threshold) | >10,000 | 1,000-10,000 | <1,000 |
| Plant TF families in top 20 hubs | >= 5 | 3-4 | <3 |
| Mean target correlation | > 0.3 | 0.15-0.3 | <0.15 |
| Genes in modules (coverage) | >80% | 50-80% | <50%|
| Known interaction recovery | >30% | 10-30% | <10% |

## Key Tools and Databases

- **PlantTFDB** (http://planttfdb.gao-lab.org/): Comprehensive plant TF annotation
  for 165 species, family classification, and binding motifs
- **PlantRegMap** (http://plantregmap.gao-lab.org/): Plant regulatory map including
  TF binding sites from ChIP-seq/DAP-seq
- **AtRegNet** (https://agris-knowledgebase.org/AtRegNet/): Curated Arabidopsis
  regulatory network
- **Plant Cistrome Database**: TF ChIP-seq data for multiple plant species
- **JASPAR Plants** (https://jaspar.genereg.net/): Plant-specific TF binding profiles

## References

- Jin et al. (2017) PlantTFDB 4.0: toward a central hub for transcription factors
  and regulatory interactions in plants. *Nucleic Acids Research*
- Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering.
  *Nature Methods*
- Huynh-Thu et al. (2010) Inferring regulatory networks from expression data using
  tree-based methods. *PLOS ONE*
- Langfelder & Horvath (2008) WGCNA: an R package for weighted correlation network
  analysis. *BMC Bioinformatics*
