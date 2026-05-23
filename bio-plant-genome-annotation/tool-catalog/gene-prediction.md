# Gene Prediction (BRAKER3)

**Goal:** Predict protein-coding genes in plant genomes using evidence-based (RNA-seq + protein hints) approach
**Best for:** Plant genome annotation with or without RNA-seq data

## Prerequisites
- BRAKER3 (includes GeneMark-ET/EP+ and AUGUSTUS)
- HISAT2 (for RNA-seq alignment)
- StringTie2 (for transcript assembly)
- ProtHint (for protein evidence)
- Repeat-masked genome (see `tool-catalog/repeat-masking.md`)

## Workflow with RNA-seq Data

### 1. Align RNA-seq Reads

```bash
# Index the masked genome
hisat2-build genome.masked.fasta genome_index

# Align paired-end reads
hisat2 -p 16 -x genome_index \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --dta | samtools sort -@ 8 -o sample_sorted.bam -
```

### 2. Assemble Transcripts

```bash
stringtie sample_sorted.bam -o sample.gtf -p 8
ls *.gtf > mergelist.txt
stringtie --merge -p 16 mergelist.txt -o merged.gtf
```

### 3. Run BRAKER3

```bash
braker.pl --genome=genome.masked.fasta \
  --bam=sample_sorted.bam \
  --hints=hints.gff \
  --species=species_name \
  --threads=16 \
  --gff3 \
  --AUGUSTUS_CONFIG_PATH=/path/to/augustus/config
```

## Workflow without RNA-seq (Protein Hints Only)

### 1. Generate Protein Hints

```bash
# ProtHint uses OrthoDB plant proteins to generate hints
prothint.py genome.masked.fasta orthodb_plant_proteins.fasta \
  --threads 16 --workdir prothint_output
```

### 2. Run BRAKER3 with Protein Hints

```bash
braker.pl --genome=genome.masked.fasta \
  --prot_seq=orthodb_plant_proteins.fasta \
  --hints=prothint_output/prothint.gff \
  --species=species_name \
  --threads=16 \
  --gff3
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| --genome | Masked genome in FASTA format |
| --bam | RNA-seq alignments (sorted BAM) |
| --prot_seq | Protein evidence (FASTA) |
| --hints | External hints file (GFF) |
| --species | AUGUSTUS species parameter set |
| --threads | CPU threads |
| --gff3 | Output GFF3 format |
| --softmasking | Use soft-masked genome |
| --fungus | Use GeneMark-ES instead of ET (not for plants) |

## Plant-Specific Considerations

- Always use `--species` with a plant-specific AUGUSTUS parameter set
- Train AUGUSTUS for your species if a close relative is not available
- For polyploid genomes: mask each subgenome separately or use subgenome-specific hints
- Plant genes have variable intron lengths — do not set strict intron size limits
- UTRs are not predicted by default in BRAKER3; use additional tools (e.g., PASA) if needed

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "AUGUSTUS_CONFIG_PATH not set" | Missing environment variable | Export path to AUGUSTUS config directory |
| "GeneMark key not found" | License issue | Install GeneMark license key in ~/.gm_key |
| "BAM file not sorted" | Unsorted alignments | Sort with `samtools sort` |
| "No genes predicted" | Hints quality too low | Check RNA-seq alignment rate; verify tissue diversity |
| "Too many short genes" | False positives from TEs | Improve repeat masking quality |
