# Bacterial and Fungal Genome Analysis Primer

## Analysis Overview

Microbial genome analysis encompasses assembly, annotation, and comparative genomics for bacteria and fungi. Unlike plant genomes, microbial genomes are smaller (bacteria: 2-12 Mb; fungi: 10-100 Mb) but require specialized tools for features like secondary metabolite gene clusters.

## Key Differences from Plant Genomes

| Feature | Bacteria | Fungi | Plants |
|---------|----------|-------|--------|
| Genome size | 2-12 Mb | 10-100 Mb | 100 Mb - 30 Gb |
| Gene structure | No introns (mostly) | Introns present | Complex introns |
| Typical genes | 2,000-8,000 | 5,000-15,000 | 20,000-50,000 |
| BGCs | Common | Common | Rare |
| Circular genome | Yes | No | No |

## Assembly Strategies

### Hybrid Assembly (Recommended)
- Combine Illumina + PacBio/ONT
- Tool: Unicycler
- Produces complete circular genomes for bacteria

### Long-read Only
- PacBio HiFi or ONT
- Tool: Flye
- Requires polishing

### Short-read Only
- Illumina paired-end
- Tool: SPAdes
- Contig-level assembly

## Annotation Components

### Gene Prediction
- **Bacteria**: Prokka (Prodigal-based)
- **Fungi**: BRAKER3 (with RNA-seq) or AUGUSTUS

### Functional Annotation
- EggNOG-mapper (COG/KEGG/GO)
- InterProScan (domains)
- SwissProt BLAST

### Specialized Annotations
- **CAZy**: Carbohydrate-active enzymes
- **antiSMASH**: Secondary metabolite clusters
- **CARD**: Antibiotic resistance genes
- **MEROPS**: Proteases
- **PHI-base**: Pathogen-host interactions

## Expected Outputs

| Output | Description |
|--------|-------------|
| Assembly FASTA | Contigs/scaffolds/chromosome |
| GFF annotation | Gene coordinates |
| Functional table | GO/KEGG/COG assignments |
| BGC report | Secondary metabolite clusters |
| Resistance report | Antibiotic resistance genes |

## Quality Metrics

| Metric | Good | Acceptable | Poor |
|--------|------|------------|------|
| N50 (bacteria) | >100 kb | >50 kb | <50 kb |
| BUSCO complete | >95% | >90% | <90% |
| Contamination | <5% | <10% | >10% |
| Completeness | >95% | >90% | <90% |

## Tools Overview

| Category | Primary Tools |
|----------|---------------|
| Assembly | Unicycler, Flye, SPAdes |
| Polishing | Pilon, Medaka, Arrow |
| Annotation | Prokka, BRAKER3 |
| BGCs | antiSMASH |
| CAZy | dbCAN2 |
| Resistance | CARD RGI |
| Comparative | Roary, FastANI, Mauve |
