# TransDecoder: Coding Region Prediction

## Tool Overview

**Tool ID**: `transdecoder`
**Category**: ORF Prediction
**Purpose**: Identify candidate coding regions in transcript sequences
**Web**: https://github.com/TransDecoder/TransDecoder
**Version**: TransDecoder v5.5+

## When to Use

- Predicting ORFs in de novo assembled transcripts
- Identifying coding potential in transcriptomes
- Preparing protein sequences for functional annotation
- Distinguishing mRNA from lncRNA

## Key Parameters

| Parameter | Default | Recommended for Plants | Description |
|-----------|---------|------------------------|-------------|
| `-t` | - | (required) | Input transcript file |
| `--min_protein_len` | 100 | **50-100** | Minimum protein length |
| `-T` | - | output_dir | Output directory |
| `--retain_pfam_hits` | - | pfam.hits | Retain ORFs with Pfam hits |
| `--retain_blastp_hits` | - | blast.hits | Retain ORFs with BLAST hits |

### Critical Parameters for Plants

#### 1. Minimum Protein Length

Default (100aa) may be too stringent for small plant proteins:

```bash
# For comprehensive ORF prediction
TransDecoder.LongOrfs -t Trinity.fasta --min_protein_len 50

# For high-confidence protein prediction
TransDecoder.LongOrfs -t Trinity.fasta --min_protein_len 100
```

#### 2. Strand-Specific Prediction

```bash
# For stranded libraries, use correct strand
TransDecoder.LongOrfs -t Trinity.fasta -S
```

## Command Examples

### Basic ORF Prediction

```bash
# Step 1: Predict long ORFs
TransDecoder.LongOrfs -t Trinity.fasta

# Step 2: Predict final ORFs
TransDecoder.Predict -t Trinity.fasta

# Output files:
# Trinity.fasta.transdecoder.pep - Predicted proteins
# Trinity.fasta.transdecoder.cds  - Predicted CDS
# Trinity.fasta.transdecoder.gff3 - ORF coordinates
# Trinity.fasta.transdecoder.bed  - BED format coordinates
```

### With Homology Evidence (Recommended)

```bash
# Step 1: Predict long ORFs
TransDecoder.LongOrfs -t Trinity.fasta

# Step 2: Search against protein database
blastp -query Trinity.fasta.transdecoder_dir/longest_orfs.pep \
  -db swissprot \
  -max_target_seqs 1 \
  -outfmt 6 \
  -evalue 1e-5 \
  -out blastp.out

# Step 3: Search against Pfam
hmmscan --cpu 8 \
  --domtblout pfam.domtblout \
  Pfam-A.hmm \
  Trinity.fasta.transdecoder_dir/longest_orfs.pep

# Step 4: Predict with homology support
TransDecoder.Predict -t Trinity.fasta \
  --retain_blastp_hits blastp.out \
  --retain_pfam_hits pfam.domtblout
```

### Plant-Specific Database Search

```bash
# Use plant-specific databases for better annotation
# Arabidopsis proteins
blastp -query longest_orfs.pep \
  -db TAIR10_pep \
  -out tair_blast.out \
  -outfmt 6 \
  -evalue 1e-5

# Plant TF database (PlantTFDB)
hmmscan --domtblout tf_domains.tbl \
  PlantTFDB.hmm \
  longest_orfs.pep
```

## Output Interpretation

### Output Files

| File | Content |
|------|---------|
| `*.transdecoder.pep` | Predicted protein sequences (FASTA) |
| `*.transdecoder.cds` | Predicted coding sequences (FASTA) |
| `*.transdecoder.gff3` | ORF coordinates in GFF3 format |
| `*.transdecoder.bed` | ORF coordinates in BED format |

### GFF3 Format

```
##gff-version 3
TRINITY_DN001_c0_g1_i1	transdecoder	gene	1	1500	.	+	.	ID=TRINITY_DN001_c0_g1_i1
TRINITY_DN001_c0_g1_i1	transdecoder	mRNA	1	1500	.	+	.	ID=TRINITY_DN001_c0_g1_i1.p1;Parent=TRINITY_DN001_c0_g1_i1
TRINITY_DN001_c0_g1_i1	transdecoder	exon	1	1500	.	+	.	ID=TRINITY_DN001_c0_g1_i1.exon1;Parent=TRINITY_DN001_c0_g1_i1.p1
TRINITY_DN001_c0_g1_i1	transdecoder	CDS	1	1500	.	+	0	ID=TRINITY_DN001_c0_g1_i1.p1.cds;Parent=TRINITY_DN001_c0_g1_i1.p1
```

### ORF Types

| Type | Description | Retained? |
|------|-------------|-----------|
| **complete** | Has start and stop codon | Yes |
| **5prime_partial** | Missing start codon | Yes |
| **3prime_partial** | Missing stop codon | Yes |
| **internal** | Missing both | Rarely |

## Plant-Specific Considerations

### 1. Alternative Splicing

Plants have extensive alternative splicing:

```bash
# One gene may have multiple isoforms with different ORFs
# Check gene_trans_map for isoform relationships

# Get gene-transcript mapping from Trinity
cat Trinity.fasta.gene_trans_map

# Compare ORFs across isoforms
```

### 2. Short Plant Proteins

Many plant proteins are short (<100aa):

```bash
# Use lower threshold for small proteins
TransDecoder.LongOrfs -t Trinity.fasta --min_protein_len 50

# But validate with homology search
```

### 3. Non-canonical Start Codons

Some plant genes use non-ATG start codons:

```bash
# TransDecoder only finds ATG-initiated ORFs
# For non-canonical starts, use additional tools
# e.g., sORF finder for small ORFs
```

## Integration with Annotation Pipeline

### Trinotate Workflow

```bash
# 1. Predict ORFs
TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta

# 2. Run BLAST against SwissProt
blastp -query Trinity.fasta.transdecoder.pep \
  -db swissprot \
  -outfmt 6 \
  -evalue 1e-5 \
  -out blastp_swissprot.out

# 3. Run Pfam search
hmmscan --cpu 8 \
  --domtblout pfam.out \
  Pfam-A.hmm \
  Trinity.fasta.transdecoder.pep

# 4. Run eggNOG-mapper
emapper.py -i Trinity.fasta.transdecoder.pep -o eggnog_out

# 5. Combine all annotations in Trinotate
```

### Separating mRNA from lncRNA

```bash
# After ORF prediction, classify transcripts

# Transcripts with ORF >= 100aa: likely mRNA
awk 'length($2) >= 300' Trinity.fasta.transdecoder.cds > mrna_candidates.cds

# Transcripts with no ORF or short ORF: potential lncRNA
# Get transcript IDs with ORFs
grep ">" Trinity.fasta.transdecoder.pep | sed 's/>//' | cut -f1 -d' ' > orf_transcripts.txt

# Subtract from all transcripts
seqkit grep -v -f orf_transcripts.txt Trinity.fasta > lncrna_candidates.fasta
```

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| No ORFs found | Wrong strand, frame | Check library strandness |
| Too many short ORFs | Low threshold | Increase `--min_protein_len` |
| Missing known proteins | Database incomplete | Use custom plant protein DB |
| Multiple ORFs per transcript | Alternative start sites | Use homology evidence to select |

## References

- TransDecoder: https://github.com/TransDecoder/TransDecoder
- Haas et al. (2013) De novo transcript sequence reconstruction from RNA-seq

## Related Tools

- **Trinity**: De novo transcriptome assembly
- **Trinotate**: Functional annotation pipeline
- **BLAST**: Protein homology search
- **HMMER/Pfam**: Domain search
- **eggNOG-mapper**: Functional annotation
