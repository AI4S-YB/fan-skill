# Plant Gene Regulatory Network: Analyst Notebook

This notebook helps you reason like an experienced plant systems biologist
when reconstructing gene regulatory networks (GRNs) from expression data.

## Before You Start: Define Your Question

GRN analysis is not one-size-fits-all. The method you choose depends
on the biological question:

- **"Which TFs regulate this pathway?"** → Supervised approach: start
  with known pathway genes, use GENIE3 or SCENIC to predict upstream TFs.
- **"What are the co-expression modules in my data?"** → Unsupervised
  approach: WGCNA, then overlay TF-target information.
- **"How does the regulatory landscape differ between conditions?"**
  → Comparative approach: build GRNs separately per condition, compare
  network topology and hub genes.
- **"Is this TF a master regulator of a biological process?"**
  → Regulon-centric approach: use SCENIC to identify the TF's regulon,
  then validate with regulon activity scores.

Clarify your question before choosing a tool. All three methods
(GENIE3, SCENIC, WGCNA) can be useful, but they answer different
questions and produce different types of output.

## GRN Inference Methods: Decision Guide

### GENIE3, SCENIC, WGCNA: When to Use Which

These three methods represent different philosophies of GRN inference.
They are not competitors -- they solve different problems.

| Dimension | GENIE3 | SCENIC | WGCNA |
|-----------|--------|--------|-------|
| Core Algorithm | Tree-based regression (Random Forest / Extra-Trees) | GENIE3 + motif enrichment (RcisTarget) + AUC scoring (AUCell) | Correlation-based co-expression + hierarchical clustering |
| Input | Expression matrix (genes × samples) | Expression matrix + TF list + motif database | Expression matrix |
| Output | Ranked TF-target predictions (weighted edge list) | Regulons (TF + direct target genes with motif support) | Co-expression modules (gene clusters) |
| Motif validation | None (purely expression-based) | Built-in (RcisTarget identifies enriched TF motifs in co-expressed gene sets) | None (purely co-expression) |
| Regulon activity quantification | Not provided | Built-in (AUCell scores per regulon per sample/cell) | Module eigengene (summary expression per module per sample) |
| Sample size requirement | n >= 15 (n >= 30 recommended) | n >= 20 (motif enrichment needs statistical power) | n >= 15 (n >= 20 for robust module detection) |
| Execution time | Medium (hours for 20K genes) | Long (days for full pipeline with many TFs) | Fast (minutes to hours) |
| Ease of biological interpretation | Moderate (long list of TF-target pairs) | High (regulons with motif support are more interpretable) | High (co-expression modules are intuitive) |
| Best for | Initial TF-target prediction, screening | Master regulator discovery, regulon definition | Co-expression module discovery, trait-module association |

**GENIE3 -- The entry point for GRN inference**:
- GENIE3 decomposes the GRN inference problem into p independent regression
  problems (one per gene). For each target gene, it uses Random Forest or
  Extra-Trees to predict its expression from the expression of all TFs.
  The importance score of each TF for predicting that target gene becomes
  the edge weight in the GRN.
- **Strength**: GENIE3 handles non-linear relationships well (tree-based
  models capture threshold effects, synergistic TF interactions). This is
  important in plants, where TF binding often depends on co-factors and
  chromatin state.
- **Weakness**: GENIE3 cannot distinguish direct (TF physically binds
  target gene promoter) from indirect (TF regulates another TF that
  regulates the target) regulation. All it sees is expression covariation.
  This is why GENIE3 edges should be treated as "regulatory hypotheses,"
  not "regulatory facts."
- **Practical tip**: GENIE3 benefits from running on TF-only expression
  as predictors (instead of all genes). By setting `regulators = tf_list`,
  you constrain the model to only consider TF-to-target edges, which
  improves both speed and interpretability.
- **Parameter**: `treeMethod = "RF"` (Random Forest) vs `"ET"` (Extra-Trees).
  Extra-Trees is faster and often performs similarly. For plant datasets
  with many highly correlated TFs (paralogs), Extra-Trees' random feature
  selection may help break ties and identify the correct paralog.

**SCENIC -- From edges to regulons with motif validation**:
- SCENIC extends GENIE3 by adding a crucial filter: motif enrichment.
  After GENIE3 generates candidate TF-target pairs (modules), RcisTarget
  scans the promoters of co-expressed target genes for enriched TF binding
  motifs. Only TF-target relationships with both expression support AND
  motif support are retained in the final regulon.
- **The motif filter is the key innovation**: it removes the majority of
  false positive edges (indirect regulation, co-expression due to
  confounding, batch effects). After motif filtering, the remaining edges
  are much more likely to represent direct regulatory relationships.
- **AUCell scores**: For each regulon, AUCell calculates an "activity score"
  per sample by ranking genes by expression and evaluating whether the
  regulon's target genes are enriched at the top of the ranking. This
  converts regulons from "static gene sets" to "dynamic activity readouts."
- **RSS (Regulon Specificity Score)**: measures how specific a regulon's
  activity is to each cell type or condition. High RSS = this regulon is
  active specifically in this context.

**WGCNA -- Co-expression modules for trait association**:
- WGCNA identifies clusters (modules) of co-expressed genes. It is NOT a
  regulatory network method per se, but modules often contain functionally
  related genes, and hub genes within modules (most connected genes) can
  be candidate regulators.
- **Best used for**: trait-module association (which co-expression module
  correlates with flowering time, yield, disease resistance?), identifying
  hub genes in modules of interest, and generating hypotheses about gene
  function ("guilt by association").
- **After WGCNA**: overlay TF-target predictions from GENIE3 or SCENIC on
  WGCNA modules. If a TF is the hub of a module AND GENIE3/SCENIC predicts
  it regulates many module members, the regulatory hypothesis is strengthened.

### Decision Flowchart

1. Do you have a list of TFs (from PlantTFDB or homology)?
   - No → Start with WGCNA (co-expression). Identify co-expression modules.
     Use homology to annotate which module genes are TFs.
   - Yes → Proceed to step 2.

2. Do you need motif support for regulatory edges (i.e., evidence beyond
   co-expression)?
   - Yes, and motif database covers your species → Use SCENIC.
   - Yes, but motif database does not cover your species → Use GENIE3
     (motif enrichment unavailable for non-model species).
   - No, initial exploratory screening is sufficient → Use GENIE3.

3. After GRN inference, do you need quantitative regulon activity per sample?
   - Yes → SCENIC (AUCell built-in).
   - No → GENIE3 edge list is sufficient for downstream analysis.

4. Do you want to associate GRN structure with a trait of interest?
   - Yes → Use WGCNA for trait-module association first, then overlay
     GENIE3/SCENIC edges on significant modules.

## TF

### Plant TF Families and Their Coverage

Understanding which TF families exist in plants, and which are well
characterized, is essential for GRN interpretation.

### MYB

The MYB superfamily is the largest TF family in plants. The R2R3-MYB
subfamily is particularly expanded (100+ members in Arabidopsis, 150+
in maize, 200+ in wheat due to polyploid expansion).

- Functions: ABA signaling, JA response, phenylpropanoid metabolism
  (flavonoids, anthocyanins), trichome/root hair development,
  secondary cell wall biosynthesis
- MYB 1R-MYB and R3-MYB subfamilies: smaller, often act as competitors
  or co-regulators with R2R3-MYBs
- Key motif: MYB binding site (MBS) variants -- type I (AACGG), type II
  (GTTAGTTA), type III (ACCAAAC)
- GRN relevance: MYBs often participate in feed-forward loops with bHLH
  and WD40 partners (the MBW complex for anthocyanin regulation)

### WRKY

WRKY TFs are defined by the conserved WRKYGQK domain and a zinc finger
motif. They bind W-box elements (TTGACC/T).

- Functions: SA-mediated defense, JA signaling crosstalk, biotic stress
  response, senescence regulation
- 70+ members in Arabidopsis, divided into 7 subgroups (I, IIa-IIe, III)
- Plant-specific family (not found in animals, only in plants and some
  protists) -- this makes them particularly interesting as plant-specific
  regulatory nodes
- GRN relevance: WRKYs often sit at the convergence points of multiple
  stress signaling pathways. A WRKY hub connecting JA and SA modules
  suggests immune signaling crosstalk.

### NAC

NAC (NAM/ATAF/CUC) TFs bind NACRS (NAC Recognition Sequence) elements.

- Functions: abiotic stress response (drought, salt, cold), leaf
  senescence, cell division, xylem differentiation, hormone signaling
- Named after founding members: NAM (No Apical Meristem), ATAF1/2,
  CUC2 (Cup-Shaped Cotyledon)
- NAC TFs are key regulators of the transition from growth to stress
  response -- they integrate environmental signals with developmental
  decisions
- GRN relevance: NAC overexpression often triggers massive transcriptomic
  reprogramming (the NAC regulon). In GRNs, NACs are often top-level
  hubs with large regulons.

### bHLH

basic helix-loop-helix TFs form dimers (homo- or hetero-) via the HLH
domain and bind E-box (CANNTG) or G-box (CACGTG) motifs via the basic region.

- Functions: JA signaling (MYC2 and related), photomorphogenesis
  (PIF family), metal homeostasis, stomatal development, flavonoid
  biosynthesis
- One of the largest TF families in plants, with extensive functional
  diversification
- The MYC2 branch of bHLH is central to JA-mediated defense and
  development. MYC2 regulons overlap with MYB and WRKY regulons in
  stress response modules.

### AP2/ERF

AP2/ERF TFs bind GCC-box (AGCCGCC, ERF subfamily) and DRE/CRT
(RCCGAC, DREB subfamily) elements.

- **DREB subfamily (CBF/DREB1)**: cold acclimation, drought, salt.
  DREB1A/CBF3 is the master regulator of cold acclimation in Arabidopsis.
  The CBF regulon includes hundreds of cold-responsive genes (COR genes).
- **ERF subfamily**: JA/ET-mediated defense. ERF1 integrates JA and ET
  signals. PDF1.2 is a classic ERF target in defense against necrotrophs.
- The DREB and ERF subfamilies represent two distinct regulatory modules
  in stress GRNs: abiotic (DREB-driven) vs biotic (ERF-driven).
- In crops (rice, wheat, maize), DREB/CBF regulons are targets for
  engineering abiotic stress tolerance. However, constitutive CBF
  overexpression often causes growth retardation (the "growth-defense
  trade-off" at the molecular level).

### bZIP

bZIP TFs bind ABRE (ABA-Responsive Element, ACGTGGC) and related motifs.

- Functions: ABA signaling core (ABI5, AREB/ABF family), seed maturation,
  osmotic stress, light signaling (HY5)
- The AREB/ABF branch is the primary mediator of ABA-dependent gene
  expression under drought stress. AREB1/ABF2, AREB2/ABF4, ABF3 are
  key ABA signaling TFs.
- bZIP TFs often form heterodimers, expanding the combinatorial
  regulatory code. GRN methods that only consider individual TFs miss
  this dimerization dimension.
- HY5 is a central integrator of light and ABA signaling, connecting
  photomorphogenesis with stress responses.

### GRAS

GRAS TFs are plant-specific and named after founding members: GAI,
RGA, SCR.

- Functions: GA signaling (DELLA proteins: GAI, RGA, RGL1-3), root
  development (SCR, SHR), nodulation (NSP1, NSP2 in legumes),
  axillary meristem formation (LAS)
- DELLA proteins are unusual: they are transcriptional repressors that
  are degraded upon GA perception. The DELLA regulatory module controls
  the GA response, which affects plant height, flowering time, and seed
  germination.
- In GRNs, DELLAs appear as "negative hubs": their degradation releases
  downstream TFs (PIFs, SPLs) from repression, activating growth programs.

### MADS-box

MADS-box TFs bind CArG-box (CC(A/T)₆GG) elements and are central to
plant developmental transitions.

- Functions: floral organ identity (ABCDE model: AP1, AP3, PI, AG, SEP),
  flowering time (FLC, SOC1, SVP), fruit development (RIN in tomato),
  seed development
- Type I MADS (SRF-like) vs Type II MADS (MEF2-like, MIKC type).
  MIKC-type MADS are the best-studied and most functionally diverse.
- MADS-box TFs often form higher-order complexes (tetramers: the "floral
  quartet model"). These protein-protein interactions are invisible to
  most GRN methods (which operate at the transcript level), but they
  determine which genes are activated or repressed.
- In crops: the MADS-box flowering time network is a major target of
  domestication and breeding. FLC and its orthologs (VRN2 in wheat/cereals)
  control vernalization requirement.

### Plant TF Databases and Their Coverage

Choosing the right TF database depends on your species:

| Database | Species Coverage | Content | URL |
|----------|-----------------|---------|-----|
| PlantTFDB v5.0 | 165 plant species | TF family classification, binding motifs (PWM), GO annotation, expression | planttfdb.gao-lab.org |
| PlantTFDB v4.0 | 120 plant species | TF family + expression atlases | planttfdb.cbi.pku.edu.cn |
| PlantRegMap | 165 species | TF binding motifs (DAP-seq, ChIP-seq), regulatory landscape | plantregmap.gao-lab.org |
| JASPAR Plants | Multiple plants | Curated TF binding motifs (PWM) | jaspar.genereg.net |
| CIS-BP | ~500 species including plants | TF binding specificities (protein binding microarray) | ccbp.ucsd.edu |
| GRASSIUS | Grasses (maize, rice, sorghum, sugarcane, Brachypodium) | Grass-specific TF families | grassius.org |
| SoyDB | Soybean | Soybean TF database | — |
| LegumeTFDB | Legumes | Legume-specific TF annotation | — |
| AGRIS | Arabidopsis | Arabidopsis TF + promoter motifs (AtcisDB) | agris-knowledgebase.org |

**Key observations for GRN analysis**:
- PlantTFDB is the primary resource for most species. If your species is
  in PlantTFDB, you can directly download the TF list for use in GENIE3
  or SCENIC (as the `regulators` or `tf_list` input).
- Motif information quality varies by species. Arabidopsis, rice, and
  maize have high-quality PWM collections. For other species, motifs are
  often transferred by orthology from Arabidopsis -- this introduces
  uncertainty (TF binding specificity can diverge between species even
  for orthologous TFs).
- For species not in PlantTFDB, you will need to infer TFs by homology
  (see "Non-Model Species TF Annotation Strategy" below).

## TC

GENIE3 tree-based

SCENIC GENIE3 +

- **RcisTarget** cis-regulatory motif
- **AUCell** regulon -- AUC
- **RSS** regulon specificity score

WGCNA TF -- WGCNA module hub

### Beyond the Three: Other GRN Methods Worth Knowing

**Partial Correlation and Graphical Models**:
- PIDC (Partial Information Decomposition) and GENENET use partial
  correlation or mutual information to distinguish direct from indirect
  edges. They model the GRN as a Gaussian Graphical Model.
- Advantage: theoretically distinguish direct vs indirect interactions
  better than GENIE3.
- Disadvantage: higher sample size requirement (n >= 50, ideally n >= 100).
  Most plant studies have insufficient samples for reliable graphical model
  estimation.

**Time-Series GRN Methods**:
- SCODE, dynGENIE3, and SINCERITIES use temporal information to infer
  directed edges (TF-A → TF-B rather than TF-A co-varies with TF-B).
- Specifically designed for time-series expression data.
- Advantage: directed edges are more mechanistically informative.
- Disadvantage: require dense time series (>= 8 time points), which limits
  applicability in many plant studies.

**Bayesian Network Methods**:
- Banjo, bnlearn in R: learn the structure of the GRN as a Directed
  Acyclic Graph (DAG) via Bayesian scoring.
- Advantage: produce a fully directed network, can incorporate prior
  knowledge (known TF-target relationships).
- Disadvantage: computationally very expensive for networks with > 100
  genes. Typically limited to focused gene sets (e.g., 50-200 candidate genes).

## Hub Gene Identification and Validation

Hub genes (highly connected nodes in the GRN) are prime candidates for
functional validation and biotechnology targets. But network degree alone
is not sufficient to claim a hub gene is biologically important.

### Hub Gene Identification

**Degree-based hubs**:
- Top 5% of genes by degree (number of edges) in the GRN
- Degree alone is a simple but noisy metric. Genes with high expression
  variance tend to have high degree in correlation-based networks
  regardless of biological relevance.

**Centrality-based hubs**:
- **Betweenness centrality**: how often a gene lies on the shortest path
  between two other genes. High betweenness = "bottleneck" genes connecting
  different network modules. These are often signaling integrators.
- **Eigenvector centrality**: a gene's importance weighted by the importance
  of its neighbors. High eigenvector centrality = connected to other highly
  connected genes (like Google PageRank for genes).
- **Closeness centrality**: average distance to all other nodes. High
  closeness = "information can flow quickly" from/to this gene.
- Combined hub score: Z-score of degree centrality + Z-score of betweenness
  centrality > 2. Genes meeting this threshold are "consensus hubs."

**Module-based hubs (WGCNA)**:
- Within each WGCNA module, the gene with the highest intramodular
  connectivity (kME) is the module hub.
- Module hubs often represent the "core function" of the module (e.g., a
  photosynthesis module hub may be a key photosynthetic enzyme or regulator).
- Module membership (kME) vs gene significance (correlation with trait):
  plotting kME vs trait correlation identifies genes that are both central
  to the module AND trait-relevant -- the strongest candidates for
  functional validation.

### Hub Gene Validation Strategy

**In silico validation**:
1. Check if the top hub TFs have been reported in the literature (PubMed
   search: "[TF name] AND [your trait/species]"). If the hub TF is known,
   it validates your network.
2. Check expression correlation: the hub TF's expression should correlate
   (r > 0.6) with its predicted targets' expression. If the TF-target
   correlations are weak, the edge may be a false positive.
3. For SCENIC regulons: check if the top hub TF's regulon is enriched
   in relevant GO terms. If the regulon targets are enriched in "defense
   response" and your experiment is a pathogen challenge, the regulon is
   biologically coherent.

**Experimental validation options** (in order of increasing effort/cost):
1. **Literature curation**: has this TF been validated (ChIP-qPCR, Y1H,
   mutant phenotype) in this or a related species? The strongest form of
   "free" validation.
2. **Public ChIP-seq / DAP-seq data**: if ChIP-seq data exists for your
   hub TF (in the same or a closely related species), check whether the
   TF binding peaks overlap with the promoters of your predicted target
   genes. This is the most direct validation method using existing data.
3. **Mutant/overexpression RNA-seq**: if an RNA-seq dataset exists for
   a mutant or overexpression line of your hub TF, check whether your
   predicted target genes are differentially expressed. This validates
   the regulatory relationship (though it cannot distinguish direct from
   indirect).
4. **Yeast One-Hybrid (Y1H)**: tests whether the TF binds to the promoter
   of a specific target gene. Medium throughput (can test dozens of
   TF-target pairs). Suitable for focused validation of top candidates.
5. **ChIP-qPCR**: tests TF binding to a specific promoter region in vivo.
   Higher confidence than Y1H because it uses the actual plant tissue.
   Low throughput (few target genes per experiment).
6. **DAP-seq**: in vitro TF binding to genomic DNA. High throughput
   (genome-wide binding sites for a single TF in one experiment). Ideal
   for validating an entire regulon at once. Requires cloning and
   expressing the TF.

### GRN Validation Approaches: Experimental Overview

**Y1H (Yeast One-Hybrid)**:
- Principle: fuse TF to GAL4 activation domain (AD). Clone target promoter
  upstream of reporter gene (HIS3 or LacZ). If TF binds promoter, reporter
  is activated, yeast grows on selective medium.
- Plant applicability: widely used for plant TFs. Works well for most TF
  families. Limitations: (a) yeast may lack plant-specific co-factors needed
  for some TF-DNA interactions; (b) chromatin context (DNA methylation,
  nucleosome positioning) absent.
- Throughput: medium (10-100 TF-target pairs tested per experiment).
- Cost: moderate.

**ChIP-qPCR**:
- Principle: cross-link TF to DNA in plant tissue, immunoprecipitate TF-DNA
  complexes using TF-specific antibody, quantify target promoter enrichment
  by qPCR.
- Plant applicability: requires either a TF-specific antibody (rare for
  most plant TFs) or a transgenic line expressing epitope-tagged TF (e.g.,
  GFP-tagged TF under native promoter). This makes ChIP-qPCR more complex
  to set up in plants than in model systems with established antibody
  resources.
- Throughput: low (1-5 target genes per experiment).
- Validation quality: high (in vivo, endogenous context).

**DAP-seq (DNA Affinity Purification sequencing)**:
- Principle: express TF in vitro, incubate with fragmented genomic DNA,
  capture TF-DNA complexes, sequence bound DNA. Genome-wide binding sites
  for a single TF.
- Plant applicability: excellent for plants. Does not require transgenic
  lines or antibodies. Has been applied to hundreds of Arabidopsis TFs
  (published as part of the cistrome project). For crop species, the
  main limitation is access to high-quality genomic DNA and the cloning
  infrastructure for expressing the TF.
- Throughput: high (genome-wide for one TF per experiment, ~1-2 weeks).
- Validation quality: high for binding specificity, but in vitro (may
  miss in vivo chromatin context effects).

## ploidy

### Polyploid Network Analysis

Polyploidy adds an extra layer of complexity to GRN inference. The
homeologous TF copies in allopolyploids or the multiple alleles in
autopolyploids create unique challenges.

- **allopolyploid** homeolog -- A/B/D genome TF
  - homeolog

- **autopolyploid** allele allelic dosage

**Allopolyploid GRN Considerations**:

1. **Homeolog-specific expression**:
   - In wheat (AABBDD), each TF gene has three homeologous copies (one
     on each subgenome: A, B, D). These copies may have diverged in
     expression pattern (tissue, stress response, developmental timing).
   - GRN inference at the gene level (summing all homeologs) collapses
     this subgenome-specific regulation into one signal, losing the
     subgenome dimension of the network.
   - If homeolog-specific expression quantification is available (requires
     subgenome-specific reads or SNP-based allele assignment), build
     separate GRNs for each subgenome to capture subgenome-specific
     regulatory logic.

2. **Homeolog cross-regulation**:
   - A TF on the A subgenome may regulate target genes on the B or D
     subgenome. This cross-subgenome regulation is invisible if subgenomes
     are analyzed separately.
   - Cross-regulating TFs often buffer the network against mutation or
     deletion in one subgenome -- a form of network robustness via
     polyploid redundancy.
   - In hexaploid wheat, cross-subgenome regulation by homoeologous TFs
     is common and contributes to the buffering capacity of the polyploid
     genome. Functional redundancy among homeologs means that single-gene
     knockouts often have no visible phenotype.

3. **Expression bias and network rewiring**:
   - After polyploidization, homeologous gene pairs may show expression
     bias (one homeolog dominates expression in certain tissues/
     conditions). This expression bias can lead to "network rewiring":
     the dominant homeolog takes over the regulatory role from the
     non-dominant one.
   - Network rewiring is a form of subfunctionalization at the regulatory
     level. It can be detected by comparing GRNs inferred from each
     subgenome's expressed TFs and targets.
   - In allopolyploid crops (wheat, cotton, oilseed rape), network
     rewiring is hypothesized to contribute to phenotypic novelty and
     adaptation to new environments.

**Autopolyploid GRN Considerations**:

- In autopolyploids (e.g., potato 4x, alfalfa 4x), each gene has
  multiple alleles rather than discrete homeologs. Allelic dosage (how
  many copies of each allele are present) can affect gene expression
  levels.
- GRN methods that use expression as input cannot distinguish between
  "TF-A allele 1 regulates Target-1" and "TF-A allele 2 regulates
  Target-2" -- they work at the gene level.
- For autopolyploids, focus analysis on gene-level GRN inference and
  acknowledge this limitation. Allele-specific GRN inference requires
  allele-specific expression quantification, which is technically
  challenging in high-ploidy autopolyploids (reads may map equally
  well to multiple alleles).

## non-model

### Non-Model Species TF Annotation Strategy

For species not covered by PlantTFDB or other TF databases, you need
a systematic approach to annotate TFs:

1. **Homology inference from Arabidopsis thaliana**:
   - Arabidopsis has the best-characterized TF repertoire of any plant.
     Use bidirectional BLAST or OrthoFinder to identify Arabidopsis
     orthologs for your species' proteins.
   - For each gene in your species: if the top Arabidopsis hit (E <
     10⁻⁵) is annotated as a TF in PlantTFDB, transfer the TF annotation.
   - Caveat: this method has false negatives (plant-specific TF families
     expanded in your species but not in Arabidopsis) and false positives
     (non-TF genes with conserved domains that match Arabidopsis TFs by
     chance).

2. **Domain-based annotation using Pfam/InterProScan**:
   - Scan all predicted proteins in your species with InterProScan.
   - Map InterPro domains to TF families using the PlantTFDB domain-family
     mapping (available in PlantTFDB's download section).
   - Advantage: identifies TFs even without close Arabidopsis orthologs.
   - Disadvantage: some domains are shared between TFs and non-TFs (e.g.,
     some kinase domains and TF domains can co-occur in the same protein).
     Manual curation may be needed for ambiguous cases.

3. **Use PlantTFDB as the primary source** (check if your species is covered):
   - PlantTFDB covers 165 species. If your species is listed, download
     the TF list directly. This is the most reliable source.

4. **Top TF hubs from known TF families as the focus**:
   - Even if complete TF annotation is challenging, focus on the major
     TF families (MYB, WRKY, NAC, bHLH, AP2/ERF, bZIP, GRAS, MADS-box)
     which are well-conserved across plants and account for the majority
     of regulatory hubs.
   - A partial but high-confidence TF list (covering these major families)
     is better than a complete but low-confidence list for GRN inference.

## GO

GO BP -- TF hub DEG

### GO Enrichment for GRN Interpretation

GO enrichment analysis is essential for interpreting GRN modules and
regulons, but it must be done correctly in the plant context:

**Background gene set selection**:
- The background should be "all expressed genes in the experiment" (genes
  that could theoretically appear in a module/regulon), NOT "all genes
  in the genome." Using the whole genome inflates significance because
  many genes are not expressed and could never appear in a module.
- Extract the background from your expression matrix: all genes with
  mean expression > threshold (e.g., count > 5 in >= 20% of samples).

**Plant GO annotation quality**:
- Arabidopsis and rice have high-quality GO annotations (experimentally
  validated + computationally inferred).
- Maize, soybean, tomato, and Medicago have moderate-quality annotations.
- Most other crops have low-quality annotations (largely inferred from
  Arabidopsis orthology, not directly validated).
- In low-quality annotation cases: (a) report this as a limitation; (b)
  use broad GO categories (GO slim) rather than specific GO terms, which
  are less affected by annotation transfer errors; (c) supplement GO with
  KEGG pathway enrichment (pathway annotations are often better curated).

**Interpreting GO enrichment for GRN hubs**:
- If a regulon's target genes are enriched in "defense response" GO:BP,
  and the hub TF is a WRKY, this is coherent and expected. Report with
  confidence.
- If a regulon is enriched in "ribosome biogenesis" and "defense response"
  simultaneously, this may indicate two distinct sub-modules within the
  regulon (the TF may regulate both growth and defense genes -- a "dual
  function" TF).
- A regulon with no GO enrichment may be: (a) false positive (low-confidence
  edges), (b) enriched in functions not captured by GO, or (c) composed
  of poorly annotated genes (common in non-model species). Interpret
  cautiously with all three possibilities in mind.

## hub

hub gene top 5% degree centrality + betweenness centrality Z-score > 2
top hub TF check expression correlation with target genes -- regulon -- in silico

### Regulon Activity Interpretation

Once you have regulons (from SCENIC) or high-confidence TF-target sets
(from GENIE3 + motif filtering), the next step is to interpret regulon
activity in biological context:

**Regulon activity as a "pathway readout"**:
- A regulon's AUCell score across samples can be treated like a
  "transcriptional pathway activity score." High AUCell = the TF's
  target genes are collectively highly expressed, suggesting the TF
  (and its regulatory program) is active in that sample.
- Compare regulon activity between conditions (e.g., drought vs control).
  If a MYB regulon has significantly higher AUCell in drought samples,
  that MYB is a drought-responsive regulator.
- Regulon activity can be more sensitive than the TF's own expression
  for detecting pathway activation, because: (a) TFs may be regulated
  post-translationally (phosphorylation, protein stability) without
  changes in their own mRNA level; (b) even a small increase in TF
  level can trigger a large change in target gene expression if the
  TF operates near its activation threshold.

**Regulon specificity across tissues/conditions**:
- Use RSS (Regulon Specificity Score) to identify regulons that are
  active only in specific contexts (e.g., flower-specific regulon,
  root-specific regulon, drought-specific regulon).
- A regulon with high RSS in one context implies that the TF's regulatory
  program is selectively deployed in that context. This is often due
  to context-specific post-translational activation of the TF or
  context-specific chromatin accessibility of its target genes.
- In crops: tissue-specific regulons are targets for precision breeding
  (modify a regulon's activity only in the target tissue, avoiding
  pleiotropic effects).

**Regulon crosstalk**:
- Calculate correlations between regulon activity scores across samples.
  If two regulons are co-active (AUCell scores correlated r > 0.8), they
  may be co-regulated or operate in the same biological program.
- Conversely, if two regulons are mutually exclusive (negatively
  correlated), they may represent antagonistic programs (e.g., growth
  vs defense, vegetative growth vs flowering).
- Regulon crosstalk analysis provides a higher-level view of the GRN,
  abstracting from individual TF-target edges to "regulatory programs."

## pitfalls

1. **TF annotation** -- co-expression false positive
2. **low sample** -- GENIE3 n >= 15 < 30 --prior
3. **batch effect** -- batch
4. **polyploid mapping** -- reads homeolog multi-mapping -- expression bias quantification

### Common Pitfalls in Plant GRN Analysis

1. **TF Annotation Gaps as a Source of False Negatives**:
   - If your TF annotation is incomplete (common in non-model species),
     true regulators may be missing from your TF list. GENIE3/SCENIC
     can only predict edges from TFs that are in the input TF list.
     Missing TFs = missing regulatory edges.
   - The most common missing TFs in non-model plants are: C2H2 zinc
     fingers, Trihelix, and TCP families -- these are less well-conserved
     and harder to annotate by homology.
   - Mitigation: in your report, note the completeness of TF annotation
     (how many TFs from each family were identified, what fraction of
     the proteome are TFs -- should be ~5-8% in plants).

2. **Co-expression as a Source of False Positives**:
   - Co-expression does not imply co-regulation. Two genes can be
     co-expressed because: (a) they share a regulatory TF, (b) they
     are in the same pathway but regulated by different TFs responding
     to the same signal, (c) they are both responsive to a common
     environmental factor without any regulatory relationship, (d)
     their expression levels are both driven by cell type composition
     in bulk tissue samples.
   - SCENIC's motif filter removes many (c) and (d) type false positives,
     but it cannot filter false positives from shared environmental
     responsiveness if the correct motif is also present by chance.
   - Mitigation: report edges with and without motif support. GENIE3
     edges without motif support = "regulatory hypotheses" (lower
     confidence). SCENIC edges with motif support = "regulatory candidates"
     (higher confidence).

3. **Low Sample Size**:
   - GENIE3: n >= 15 is the absolute minimum, but standard errors of
     edge weights are large. n >= 30 is recommended for stable edge
     rankings (top edges more likely to replicate in independent data).
   - n < 30: the top 10% of edges (by weight) are relatively stable;
     edges below the top 30% may shuffle completely with a different
     sample. Do not overinterpret moderate-weight edges.
   - If n < 15, focus on prior-supported analysis: pre-select a set of
     candidate TFs based on literature/functional annotation, and only
     evaluate these TFs' predicted targets (rather than all TFs).
     Report this as "candidate TF-target analysis" rather than "GRN
     inference."
   - Impact on SCENIC: with low n, the motif enrichment step (RcisTarget)
     loses power because co-expression modules are less robust. Fewer
     regulons pass the motif filter. Report the fraction of GENIE3 modules
     that passed RcisTarget filtering.

4. **Batch Effects Creating Spurious Edges**:
   - If samples from batch 1 and batch 2 have systematic expression
     differences, genes that are "batch-responsive" will cluster together
     and appear co-expressed. GENIE3 will infer regulatory edges among
     them -- a "batch regulon" that is purely technical.
   - Detection: after batch correction, rebuild the GRN. If the top hub
     genes and the top edges change substantially (> 50% of top edges
     are different), batch effects were driving the original network.
   - If batch correction cannot fully remove batch effects (e.g., batch
     is confounded with treatment), report this limitation and restrict
     network interpretation to edges that persist after batch correction.

5. **Polyploid Mapping Issues**:
   - In RNA-seq of polyploid species, reads from homeologous genes may
     multi-map (align equally well to multiple subgenomes). This inflates
     expression estimates for some homeologs and deflates them for others.
   - In GRN inference, if homeolog-A is overestimated (due to cross-mapping
     of homeolog-B reads), it may appear to be a hub with many connections
     -- these connections are artefacts of the inflated expression variance.
   - Mitigation: use only uniquely mapping reads for expression
     quantification in polyploid species. If uniquely mapping rate is low
     (< 50%), homeolog-level expression estimates are unreliable. Use
     gene-family-level rather than homeolog-level quantification.

### Indirect Regulation vs Direct Binding

A fundamental limitation of expression-based GRN inference: it cannot
distinguish a TF that directly binds and activates a target gene from
a TF that activates an intermediate TF which in turn activates the
target. Both scenarios produce the same expression covariation pattern.

**How indirect regulation manifests in GRNs**:
- In a three-gene chain: TF-A → TF-B → Target-C, GENIE3 may assign high
  weight to both edges (A→B and A→C). The A→C edge is indirect but
  appears strong because A's expression correlates with C's expression
  (through B). A→C is a "transitive edge" (false direct edge).
- Transitive edges inflate the degree of upstream TFs, making them appear
  as massive hubs (regulating hundreds of genes) when many of those edges
  are indirect.

**Mitigation strategies**:
1. Motif filtering (SCENIC): if the A→C edge lacks a TF-A binding motif
   in C's promoter, SCENIC will remove it even if the expression correlation
   is strong. This is the most effective filter for indirect edges.
2. Partial correlation / graphical models: these methods condition on
   intermediate genes. If conditioning on B's expression removes the
   correlation between A and C, the A→C edge is likely indirect.
3. Time-series data: in a time course, A changes first, then B, then C.
   If A→C is direct, the time lag A→C should be shorter than A→B→C.
   Time-series GRN methods exploit this.

### Motif Enrichment False Positives

Motif enrichment (RcisTarget in SCENIC) is a powerful filter, but it
has its own false positive risks in plants:

**Short, degenerate motifs**:
- Some TF binding motifs are short (6-8 bp) and degenerate (multiple
  possible bases at some positions). Such motifs occur by chance
  thousands of times in the genome.
- If RcisTarget finds an enrichment of a short degenerate motif in a
  co-expressed gene set, it may be a statistical artefact (the motif is
  so common that any random gene set of that size would be "enriched").
- Check the motif's information content (IC). Motifs with IC < 8 bits
  are low-complexity and prone to false positive enrichment.

**GC content bias**:
- Plant promoters can have extreme GC content (from < 30% to > 60%).
- A co-expressed gene set may share similar GC content (e.g., highly
  expressed genes tend to have higher GC content in some species).
- RcisTarget's background model assumes gene sets are random with respect
  to sequence composition. If the co-expressed gene set has biased GC
  content, motifs matching that GC bias will appear enriched -- even if
  they have no biological relevance.
- Mitigation: use a GC-matched background (genes with similar GC content
  as the co-expressed set) for motif enrichment testing.

**Conserved non-coding sequences mistaken for functional motifs**:
- Conserved non-coding sequences (CNSs) in plant promoters may be conserved
  for reasons other than TF binding (e.g., structural DNA features,
  nucleosome positioning). RcisTarget may identify a motif within a CNS
  and call it enriched, but the motif may not be the functional driver
  of the CNS's conservation.
- Cross-reference enriched motifs with known TF binding data (DAP-seq,
  ChIP-seq). If a motif is enriched in your gene set AND known to be
  bound by a TF in your species (from public data), the functional
  relevance is much stronger.
