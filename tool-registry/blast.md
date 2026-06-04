# BLAST: Sequence Similarity Search

## Tool Overview

**Tool ID**: `blast`
**Category**: Sequence Similarity Search
**Purpose**: Identify homologous sequences through local alignment
**Source**: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
**Version**: BLAST+ 2.x (recommended)

## When to Use

- Gene family member identification and validation
- Homology search against reference databases
- Cross-species sequence comparison
- Synteny analysis input generation (for MCScanX)
- Annotation transfer between species

## BLAST Programs

| Program | Query | Database | Use Case |
|---------|-------|----------|----------|
| **blastp** | Protein | Protein | Gene family identification, homolog search |
| **blastn** | Nucleotide | Nucleotide | DNA sequence search |
| **blastx** | Nucleotide (translated) | Protein | Gene prediction from DNA |
| **tblastn** | Protein | Nucleotide (translated) | Find genes in unannotated genomes |
| **tblastx** | Nucleotide (translated) | Nucleotide (translated) | Sensitive DNA-DNA comparison |

## Key Parameters

| Parameter | Default | Recommended | Description |
|-----------|---------|-------------|-------------|
| `-evalue` | 10 | **1e-5 to 1e-10** | E-value threshold |
| `-outfmt` | 0 | **6** | Tabular output format |
| `-max_target_seqs` | 500 | **5-20** | Max hits per query |
| `-qcov_hsp_perc` | - | **50-70** | Query coverage threshold |
| `-num_threads` | 1 | **8-16** | Number of CPU threads |

### Critical Parameters for Gene Family Analysis

#### 1. E-value Threshold

```bash
# Gene family identification: use stringent e-value
blastp -query query.fasta -db proteins.fasta -evalue 1e-10 -outfmt 6

# Sensitive search (remote homologs): use lenient e-value
blastp -query query.fasta -db proteins.fasta -evalue 1e-5 -outfmt 6
```

#### 2. Coverage Filtering

```bash
# Filter by query coverage and subject coverage
# Column 3: % identity, Column 4: alignment length
# Column 13: query length, Column 14: subject length

awk '{
  qcov = $4/$13*100;
  scov = $4/$14*100;
  if (qcov >= 50 && scov >= 50 && $3 >= 30) print
}' blast_results.out > filtered_results.out
```

## Command Examples

### Gene Family Identification (Primary)

```bash
# Step 1: Create BLAST database
makeblastdb -in genome_proteins.fasta -dbtype prot -out proteins_db

# Step 2: Run BLASTP with known family members
blastp -query known_family_members.fasta \
  -db proteins_db \
  -evalue 1e-10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -max_target_seqs 1000 \
  -num_threads 8 \
  -out blast_results.out

# Step 3: Filter by identity and coverage
awk '$3 >= 30 && ($4/$13) >= 0.5 && ($4/$14) >= 0.5' blast_results.out > filtered_hits.out

# Step 4: Extract hit sequences
cut -f2 filtered_hits.out | sort -u | \
  seqkit grep -f - genome_proteins.fasta > family_candidates.fasta
```

### Iterative Search Strategy (HMM-like)

```bash
# Round 1: Initial BLAST with seed sequences
blastp -query seed_sequences.fasta \
  -db proteins_db \
  -evalue 1e-5 \
  -outfmt 6 \
  -out round1.hits

# Round 2: Use hits as new queries
cut -f2 round1.hits | sort -u > hit_ids.txt
seqkit grep -f hit_ids.txt proteins_db.fasta > round2_queries.fasta

blastp -query round2_queries.fasta \
  -db proteins_db \
  -evalue 1e-5 \
  -outfmt 6 \
  -out round2.hits

# Combine and deduplicate
cat round1.hits round2.hits | cut -f2 | sort -u > all_family_members.txt
```

### MCScanX Input Generation

```bash
# All-vs-all BLAST for synteny analysis
blastp -query proteins.fasta \
  -db proteins.fasta \
  -evalue 1e-5 \
  -outfmt 6 \
  -max_target_seqs 5 \
  -num_threads 16 \
  -out gene_location.blast

# This output is used directly by MCScanX
```

### Cross-Species Homolog Search

```bash
# Find orthologs in related species
blastp -query species1_genes.fasta \
  -db species2_proteins.fasta \
  -evalue 1e-10 \
  -outfmt 6 \
  -max_target_seqs 10 \
  -out cross_species_blast.out

# Find best hits (reciprocal best hit strategy)
awk '!seen[$1]++' cross_species_blast.out > best_hits.out
```

## Output Format 6 Columns

| Column | Name | Description |
|--------|------|-------------|
| 1 | qseqid | Query sequence ID |
| 2 | sseqid | Subject (hit) sequence ID |
| 3 | pident | Percentage identity |
| 4 | length | Alignment length |
| 5 | mismatch | Number of mismatches |
| 6 | gapopen | Number of gap opens |
| 7 | qstart | Query alignment start |
| 8 | qend | Query alignment end |
| 9 | sstart | Subject alignment start |
| 10 | send | Subject alignment end |
| 11 | evalue | E-value |
| 12 | bitscore | Bit score |

### Extended Format (Recommended)

```bash
# Include query and subject lengths for coverage calculation
blastp -query query.fasta -db db.fasta \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -out blast_results.out
```

## Plant-Specific Considerations

### 1. Multi-Gene Family Search

Plant gene families often have many paralogs:

```bash
# Use higher max_target_seqs for large families
blastp -query query.fasta \
  -db plant_proteins.fasta \
  -max_target_seqs 500 \
  -evalue 1e-5 \
  -outfmt 6

# Then cluster results to remove redundancy
cd-hit -i family_candidates.fasta -o clustered.fasta -c 0.9
```

### 2. Polyploid Species

```bash
# Search each subgenome separately first
blastp -query query.fasta -db subgenome_A.fasta -out A_hits.out
blastp -query query.fasta -db subgenome_D.fasta -out D_hits.out

# Then combine and compare homeologs
cat A_hits.out D_hits.out > all_hits.out
```

### 3. Distant Homolog Detection

For detecting remote homologs in divergent plant species:

```bash
# Use lower e-value threshold and PSI-BLAST for distant homologs
psiblast -query query.fasta \
  -db proteins.fasta \
  -evalue 1e-3 \
  -num_iterations 3 \
  -outfmt 6 \
  -out psiblast_results.out
```

## Integration with Gene Family Workflow

```bash
# Complete BLAST validation workflow

# 1. HMM search (primary method)
hmmsearch --tblout hmm_hits.tbl SBP.hmm proteins.fasta

# 2. BLAST validation (secondary method)
blastp -query known_SPL.fasta \
  -db proteins.fasta \
  -evalue 1e-10 \
  -outfmt 6 \
  -out blast_hits.out

# 3. Combine results (intersection for high confidence)
cut -f1 hmm_hits.tbl | grep -v "^#" | sort -u > hmm_ids.txt
cut -f2 blast_hits.out | sort -u > blast_ids.txt
comm -12 hmm_ids.txt blast_ids.txt > high_confidence_ids.txt

# 4. Union for comprehensive coverage
cat hmm_ids.txt blast_ids.txt | sort -u > all_candidate_ids.txt
```

## Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| Too many hits | E-value too lenient | Use `-evalue 1e-10` or stricter |
| Missing known homologs | Database incomplete | Check BUSCO completeness |
| Low identity hits | Distant homologs | Use PSI-BLAST or HMM methods |
| Memory error | Large database | Use `-remote` for NCBI or split queries |
| Slow search | Too many threads | Reduce to 8-16 threads |

## References

- Altschul et al. (1990) Basic local alignment search tool
- Camacho et al. (2009) BLAST+: architecture and applications
- NCBI BLAST+: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

## Related Tools

- **HMMER**: More sensitive for remote homologs
- **DIAMOND**: Faster alternative to BLAST for large datasets
- **MMseqs2**: Ultra-fast sequence search
- **CD-HIT**: Cluster sequences by similarity
