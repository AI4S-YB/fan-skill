# Bacterial and Fungal Genome Analysis Notebook

## Overview

This notebook provides step-by-step guidance for bacterial and fungal genome assembly, annotation, and comparative genomics.

---

## Phase 1: Quality Control

### 1.1 Assess Read Quality

```bash
# FastQC for all samples
fastqc *.fastq.gz -o fastqc_results/

# MultiQC summary
multiqc fastqc_results/ -o multiqc_report/
```

### 1.2 Adapter Trimming

```bash
# Using fastp
fastp -i R1.fastq.gz -I R2.fastq.gz \
    -o R1.trimmed.fastq.gz -O R2.trimmed.fastq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --html fastp.html
```

### 1.3 Long-read Quality (if applicable)

```bash
# NanoPlot for ONT data
NanoPlot --fastq ont_reads.fastq.gz \
    -o nanoplot_results/

# Filter long reads by length and quality
Filtlong --min_length 1000 \
    --mean_q 10 \
    --length_percentile 90 \
    ont_reads.fastq.gz > filtered_reads.fastq
```

---

## Phase 2: Genome Assembly

### 2.1 Hybrid Assembly with Unicycler

```bash
# Recommended for Illumina + long-read data
unicycler -1 R1.trimmed.fastq.gz -2 R2.trimmed.fastq.gz \
    -l filtered_reads.fastq \
    -o unicycler_output/ \
    -t 32 \
    --mode bold

# Output: assembly.fasta
```

### 2.2 Long-read Assembly with Flye

```bash
# For long-read only data
flye --nano-raw ont_reads.fastq.gz \
    --genome-size 5m \
    --out-dir flye_output/ \
    --threads 32 \
    --iterations 2

# For PacBio HiFi
flye --pacbio-hifi hifi_reads.fastq.gz \
    --genome-size 5m \
    --out-dir flye_output/ \
    --threads 32
```

### 2.3 Short-read Assembly with SPAdes

```bash
# For Illumina only data
spades.py -1 R1.trimmed.fastq.gz -2 R2.trimmed.fastq.gz \
    -o spades_output/ \
    --careful \
    --cov-cutoff auto \
    -t 32 -m 128
```

---

## Phase 3: Polishing

### 3.1 Medaka (for ONT data)

```bash
# ONT-specific polishing
medaka_consensus -i ont_reads.fastq.gz \
    -d assembly.fasta \
    -o medaka_output \
    -t 16 \
    -m r941_min_high_g360

# Output: medaka_output/consensus.fasta
```

### 3.2 Pilon (using Illumina data)

```bash
# Map Illumina reads to assembly
bwa index polished_assembly.fasta
bwa mem -t 16 polished_assembly.fasta R1.fastq.gz R2.fastq.gz | \
    samtools sort -o aligned.bam
samtools index aligned.bam

# Run Pilon
pilon --genome polished_assembly.fasta \
    --bam aligned.bam \
    --output pilon_polished \
    --threads 16 \
    --fix all

# Typically 2-3 rounds needed
```

### 3.3 Arrow/Gcpp (for PacBio data)

```bash
# Using gcpp (GenomicConsensus)
gcpp --verbose --algorithm=arrow \
    --reference assembly.fasta \
    -o arrow_polished.fasta \
    subreads.bam
```

---

## Phase 4: Assembly Quality Assessment

### 4.1 QUAST

```bash
quast.py assembly.fasta \
    -o quast_output/ \
    --min-contig 500 \
    --gene-finding \
    -t 16
```

### 4.2 BUSCO

```bash
# For bacteria
busco -i assembly.fasta \
    -l bacteria_odb10 \
    -o busco_output \
    -m genome \
    -c 16

# For fungi
busco -i assembly.fasta \
    -l fungi_odb10 \
    -o busco_output \
    -m genome \
    -c 16
```

### 4.3 CheckM (for bacteria)

```bash
checkm lineage_wf assembly.fasta checkm_output/ \
    -t 16 --pplacer_threads 8

# Output includes completeness and contamination estimates
```

---

## Phase 5: Gene Annotation

### 5.1 Prokka (Bacteria)

```bash
# Standard bacterial annotation
prokka --outdir prokka_output \
    --prefix strain_name \
    --genus Escherichia \
    --species coli \
    --strain K12 \
    --kingdom Bacteria \
    --usegenus \
    --cpus 16 \
    assembly.fasta

# Output includes: .gff, .faa, .ffn, .gbk
```

### 5.2 BRAKER3 (Fungi with RNA-seq)

```bash
# Fungal annotation with RNA-seq evidence
braker.pl --genome=assembly.fasta \
    --bam=rna_seq_aligned.bam \
    --species=fungus_name \
    --softmasking \
    --cores=32 \
    --etpmode

# Use fungi BUSCO lineage
```

### 5.3 AUGUSTUS (Fungi without RNA-seq)

```bash
# Use pre-trained fungal model
augustus --species=aspergillus_nidulans \
    --gff3=on \
    assembly.fasta > genes.gff3

# Available fungal species models:
# aspergillus_nidulans, candida_albicans, saccharomyces_cerevisiae, etc.
```

---

## Phase 6: Functional Annotation

### 6.1 EggNOG-mapper

```bash
# Run EggNOG-mapper on protein sequences
emapper.py -i proteins.faa \
    -o eggnog_output \
    --data_dir /path/to/eggnog_db \
    --tax_scope 2  # Bacteria
    --go_evidence \
    --cpu 16

# For fungi, use --tax_scope 4751
```

### 6.2 InterProScan

```bash
interproscan.sh -i proteins.faa \
    -o interpro_output.tsv \
    -f TSV,XML,GFF3 \
    -dp \
    -goterms \
    -pa \
    -cpu 16
```

---

## Phase 7: CMPP Specialized Annotations

### 7.1 CAZy (Carbohydrate-active enzymes)

```bash
# Using dbCAN2
run_dbcan proteins.faa protein \
    --out_dir cazy_output \
    --db_dir /path/to/dbcan_db

# Outputs: CAZy families (GH, GT, PL, CE, CBM, AA)
```

### 7.2 MEROPS (Proteases)

```bash
# BLAST against MEROPS database
blastp -query proteins.faa \
    -db merops \
    -evalue 1e-5 \
    -outfmt 6 \
    -out merops_hits.tsv

# Filter by identity >= 30%
awk '$3 >= 30' merops_hits.tsv > filtered_merops.tsv
```

### 7.3 P450 Enzymes

```bash
# HMM search against P450 database
hmmscan --cpu 16 \
    --tblout p450_hits.tsv \
    p450.hmm proteins.faa
```

### 7.4 PHI-base (Pathogen-host)

```bash
# BLAST against PHI-base
blastp -query proteins.faa \
    -db phi_base \
    -evalue 1e-5 \
    -max_target_seqs 5 \
    -outfmt "6 qseqid sseqid pident length evalue stitle" \
    -out phi_hits.tsv
```

---

## Phase 8: Secondary Metabolite Prediction

### 8.1 antiSMASH (Bacteria)

```bash
# Bacterial BGC prediction
antismash assembly.fasta \
    --output-dir antismash_output \
    --taxon bacteria \
    --genefinding-tool prodigal \
    --cb-knownclusters \
    --cb-subclusters \
    --asf \
    --smcog-trees \
    --cpus 16

# Output includes: BGC types, known cluster matches, core enzyme predictions
```

### 8.2 antiSMASH (Fungi)

```bash
# Fungal BGC prediction
antismash assembly.fasta \
    --output-dir antismash_output \
    --taxon fungi \
    --genefinding-tool none \  # Use existing GFF
    --genefinding-gff3 genes.gff3 \
    --cb-knownclusters \
    --cpus 16
```

---

## Phase 9: Antibiotic Resistance Genes

### 9.1 CARD RGI

```bash
# Run CARD resistance gene identifier
rgi main --input_sequence assembly.fasta \
    --output_file rgi_output \
    --input_type contig \
    --clean \
    --num_threads 16

# Output: resistance genes with phenotype predictions
```

### 9.2 ResFinder

```bash
# Alternative resistance finder
resfinder.py -i assembly.fasta \
    -o resfinder_output \
    -s ecoli \
    -t 0.9  # 90% identity threshold
```

---

## Phase 10: Comparative Genomics

### 10.1 Pan-genome with Roary

```bash
# Create input directory with GFF files
mkdir gff_files/
# Copy all .gff files to gff_files/

# Run Roary
roary -e -n -p 16 \
    -f roary_output \
    -i 95 \  # BLAST identity threshold
    -cd 99 \  # Core genes in 99% of isolates
    gff_files/*.gff

# Output: core/pan gene lists, gene presence matrix, phylogeny
```

### 10.2 ANI Calculation

```bash
# FastANI for genome similarity
fastani -q genome1.fasta \
    -r genome2.fasta \
    -o ani_result.txt

# Pairwise ANI for multiple genomes
fastANI.py --query_list genomes.txt \
    --ref_list genomes.txt \
    -o ani_matrix.tsv \
    -t 16
```

### 10.3 Core-genome Phylogeny

```bash
# Using IQ-TREE on core alignment from Roary
iqtree2 -s roary_output/core_gene_alignment.aln \
    -m MFP \
    -bb 1000 \
    -nt 16

# Output: phylogeny tree with bootstrap support
```

### 10.4 Synteny Analysis

```bash
# Using Mauve
progressiveMauve --output=alignment.xmfa \
    genome1.fasta genome2.fasta

# Or using MUMmer
nucmer --maxmatch -p nucmer_output \
    reference.fasta query.fasta
mummerplot nucmer_output.delta --png
```

---

## Phase 11: Visualization

### 11.1 Circular Genome Plot

```bash
# Using CGView
cgview -i genes.gff \
    -s assembly.fasta \
    -o genome_map.png \
    --size 1000

# Or using DNAPlotter (GUI)
# Or using Proksee (web-based)
```

### 11.2 BGC Visualization

```bash
# antiSMASH outputs HTML visualization
# Open antismash_output/index.html in browser
```

---

## FAQ

### Q1: Assembly is fragmented (many contigs)

**Solutions:**
- Use hybrid assembly if long reads available
- Increase sequencing depth
- Check for contamination
- Use different k-mer sizes in SPAdes

### Q2: BUSCO completeness is low

**Solutions:**
- Check for contamination
- Verify correct lineage selected
- Check assembly for gaps
- May need more data

### Q3: antiSMASH finds no clusters

**Possible causes:**
- Strain truly lacks BGCs
- Assembly incomplete
- Gene prediction failed
- Try different gene finding tool

### Q4: Plasmid assembly issues

**Solutions:**
- Use plasmidSPAdes for short reads
- Check coverage (plasmids may have high copy)
- May need manual curation for circularization
