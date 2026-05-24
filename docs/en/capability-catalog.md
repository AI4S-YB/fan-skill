# Fan-Skill Analysis Capability Catalog

30 analysis entries organized by what you want to achieve.

## Find Genes Controlling Traits

| Entry | What it does | Typical input |
|------|------|------|
| `gwas` | Genome-wide association — scan for trait-associated loci | VCF + phenotype |
| `qtl-mapping` | Linkage-based QTL for biparental populations | Genetic map + phenotype |
| `candidate-gene-association` | Association within pre-selected gene regions | VCF + gene list + phenotype |
| `eqtl` | Expression QTL — variants associated with expression | Genotype + expression matrix |
| `population` | Population structure — PCA, ADMIXTURE, Fst, phylogeny | VCF/PLINK |

## Understand Gene Function

| Entry | What it does | Typical input |
|------|------|------|
| `rnaseq` | Differential expression | Count matrix + metadata |
| `time-series` | Developmental/time-course patterns | Expression + time points |
| `grn` | Gene regulatory network inference | Expression + TF list |
| `small-rna` | miRNA prediction, target prediction | sRNA-seq reads + reference |
| `multi-omics` | Multi-omics integration | Multiple omics datasets |

## Predict Breeding Value

| Entry | What it does | Typical input |
|------|------|------|
| `genomic-selection` | Predict breeding values | Genotype + phenotype |
| `hybrid-prediction` | Hybrid performance, GCA/SCA | Parental genotypes + hybrid phenotypes |
| `phenotype` | Heritability, BLUP/BLUE, MET summary | Phenotype + trial design |
| `enviromics` | Environmental + genomic integration | Climate + genotype + phenotype |

## Data Processing

| Entry | What it does | Typical input |
|------|------|------|
| `variant-calling` | SNP/InDel discovery | FASTQ + reference |
| `genotype-imputation` | Fill missing genotypes | VCF + reference panel |
| `genome-assembly` | De novo assembly (HiFi, ONT, hybrid) | Sequencing reads |
| `genome-annotation` | Gene prediction, functional annotation | Assembly + RNA-seq evidence |

## Epigenomics & Chromatin

| Entry | What it does | Typical input |
|------|------|------|
| `chipseq` | ChIP-seq peak calling, differential binding | ChIP + input FASTQ |
| `atacseq` | Open chromatin identification | ATAC-seq FASTQ |
| `methylation` | DNA methylation (WGBS/RRBS) | Bisulfite-seq reads |

## Microbiome & Metagenomics

| Entry | What it does | Typical input |
|------|------|------|
| `amplicon` | 16S/ITS amplicon analysis | Amplicon FASTQ |
| `metagenomics` | Shotgun assembly, binning, MAGs | Metagenomic reads |

## Metabolomics & Proteomics

| Entry | What it does | Typical input |
|------|------|------|
| `metabolomics` | LC-MS/GC-MS, differential metabolites | Mass spectrometry data |
| `proteomics` | DDA/DIA quantification, differential proteins, PPI | Mass spectrometry data |

## Evolution & Comparative Genomics

| Entry | What it does | Typical input |
|------|------|------|
| `comparative` | Synteny, Ks, gene families, selection pressure | Genome + annotation |
| `pan-genome` | Core/variable genome, PAV detection | Multiple assemblies |

## Genome Editing & Markers

| Entry | What it does | Typical input |
|------|------|------|
| `crispr` | sgRNA design, off-target prediction | Target gene sequence |
| `marker` | KASP/InDel/SSR development, parental recommendation | GWAS peak / QTL interval |

## Visualization

| Entry | What it does | Typical input |
|------|------|------|
| `visualization` | Publication-quality figures | Analysis output |

## Plant-Specific Capabilities

- **Polyploid-aware**: auto-detection, subgenome-specific analysis
- **Self/cross-pollinated**: different GWAS and GS strategies
- **12 crop species cheatsheet**: reference genomes, LD decay, breeding systems
- **Non-model species strategy**: Mercator annotation, cross-species inference
- **Multi-environment trial (MET)**: standard in plant breeding
