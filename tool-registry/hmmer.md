# HMMER: Profile Hidden Markov Model Search

## Tool Overview

**Tool ID**: `hmmer`
**Category**: Sequence Analysis / Homology Search
**Purpose**: Sensitive protein sequence analysis using profile hidden Markov models
**Source**: http://hmmer.org/
**Version**: HMMER 3.3+ (recommended)

## When to Use

- Gene family member identification (primary method)
- Protein domain detection
- Remote homolog discovery
- More sensitive than BLAST for divergent sequences

## Key Programs

| Program | Purpose |
|---------|---------|
| **hmmsearch** | Search HMM profile against sequence database |
| **hmmscan** | Search sequence against HMM database (Pfam) |
| **hmmbuild** | Build HMM profile from multiple alignment |
| **hmmalign** | Align sequences to HMM profile |
| **jackhmmer** | Iterative search (like PSI-BLAST) |

## Key Parameters

### hmmsearch

| Parameter | Default | Recommended | Description |
|-----------|---------|-------------|-------------|
| `-E` | 10 | **1e-5** | E-value threshold for full sequence |
| `--domE` | 10 | **1e-5** | E-value threshold for domains |
| `--tblout` | - | Required | Table output file |
| `--domtblout` | - | Optional | Domain-level output |
| `--cpu` | 2 | **8** | Number of threads |
| `--cut_tc` | - | Use with Pfam | Use trusted cutoff from HMM |

### Critical Parameters for Gene Family Analysis

#### 1. E-value Threshold

```bash
# Standard gene family search
hmmsearch -E 1e-5 --tblout results.tbl domain.hmm proteins.fasta

# Sensitive search for remote homologs
hmmsearch -E 1e-3 --tblout results.tbl domain.hmm proteins.fasta
```

#### 2. Coverage/Domain Check

```bash
# Check both sequence and domain E-values
hmmsearch -E 1e-5 --domE 1e-5 \
  --tblout seq_hits.tbl \
  --domtblout domain_hits.tbl \
  domain.hmm proteins.fasta
```

## Command Examples

### Gene Family Identification (Primary Workflow)

```bash
# Step 1: Get Pfam HMM profile
# Option A: Download from Pfam
wget https://pfam.xfam.org/family/PF03110/hmm -O SBP.hmm

# Option B: From local Pfam installation
hmmfetch Pfam-A.hmm PF03110 > SBP.hmm

# Step 2: Search against genome proteins
hmmsearch --cpu 8 \
  -E 1e-5 \
  --tblout hmm_results.tbl \
  SBP.hmm \
  genome_proteins.fasta

# Step 3: Extract hit IDs
grep -v "^#" hmm_results.tbl | awk '{print $1}' | sort -u > family_ids.txt

# Step 4: Extract sequences
seqkit grep -f family_ids.txt genome_proteins.fasta > family_proteins.fasta
```

### Building Custom HMM Profile

```bash
# Step 1: Align known family sequences
mafft --auto known_family_members.fasta > aligned.fasta

# Step 2: Build HMM profile
hmmbuild custom_family.hmm aligned.fasta

# Step 3: Search with custom profile
hmmsearch -E 1e-5 --tblout results.tbl \
  custom_family.hmm proteins.fasta

# Step 4: Optionally calibrate for better E-values
hmmcalibrate custom_family.hmm  # Takes time but improves sensitivity
```

### Iterative Search (jackhmmer)

```bash
# Iterative search for comprehensive coverage
jackhmmer -N 5 \
  -E 1e-5 \
  --tblout jackhmmer_results.tbl \
  seed_sequence.fasta \
  proteins.fasta

# -N 5: 5 iterations maximum
# Useful for finding remote homologs
```

### Pfam Domain Annotation

```bash
# Scan proteins against all Pfam domains
hmmscan --cpu 8 \
  --domtblout pfam_domains.tbl \
  --cut_tc \
  Pfam-A.hmm \
  proteins.fasta

# --cut_tc: Use trusted cutoff for each Pfam domain
```

## Output Format

### Table Output (--tblout)

```
# target name        accession query name accession E-value  score bias ...
Gene001              -        SBP        PF03110   1.2e-45  150.2 0.0
Gene002              -        SBP        PF03110   3.4e-38  128.5 0.0
```

### Domain Output (--domtblout)

```
# target name  accession tlen query name accession qlen E-value score ...
Gene001        -         500  SBP        PF03110   78   1.2e-45 150.2 ...
  == domain 1 of 1: 120-398 (target) matches 5-72 (query)
```

## Integration with BLAST Validation

```bash
# Combined HMM + BLAST approach for high confidence

# 1. HMM search (primary)
hmmsearch -E 1e-5 --tblout hmm_hits.tbl SBP.hmm proteins.fasta
cut -f1 hmm_hits.tbl | grep -v "^#" | sort -u > hmm_ids.txt

# 2. BLAST search (validation)
blastp -query known_SPL.fasta -db proteins.fasta \
  -evalue 1e-10 -outfmt 6 -out blast_hits.out
cut -f2 blast_hits.out | sort -u > blast_ids.txt

# 3. Intersection (high confidence)
comm -12 hmm_ids.txt blast_ids.txt > high_confidence_ids.txt

# 4. Union (comprehensive)
cat hmm_ids.txt blast_ids.txt | sort -u > all_candidate_ids.txt
```

## Plant-Specific Considerations

### 1. Multi-Domain Proteins

Plant proteins often have multiple domains:

```bash
# Use --domtblout to see all domain hits
hmmsearch --domtblout all_domains.tbl \
  Pfam-A.hmm plant_proteins.fasta

# Filter for multi-domain proteins
awk '{print $1}' all_domains.tbl | sort | uniq -c | \
  awk '$1 > 1' > multi_domain_proteins.txt
```

### 2. Large Gene Families

```bash
# For families with 100+ members, use efficient parameters
hmmsearch --cpu 16 \
  --max        * Report all hits, not just top *
  -E 1e-5 \
  --tblout results.tbl \
  domain.hmm large_proteome.fasta
```

### 3. Polyploid Genomes

```bash
# Search each subgenome separately
hmmsearch -E 1e-5 --tblout A_hits.tbl SBP.hmm subgenome_A.fasta
hmmsearch -E 1e-5 --tblout D_hits.tbl SBP.hmm subgenome_D.fasta

# Combine results
cat A_hits.tbl D_hits.tbl | grep -v "^#" > all_hits.tbl
```

## Performance Tips

| Scenario | Recommendation |
|----------|----------------|
| Small database (<10K seqs) | Standard hmmsearch |
| Large database (>100K seqs) | Use `--cpu 16`, consider `--noali` |
| Multiple HMMs | Use hmmpress + hmmscan |
| Sensitive search | Use jackhmmer (iterative) |
| Speed priority | Use `--F1 0.02 --F2 0.001` |

## Common Issues

| Problem | Cause | Solution |
|---------|-------|----------|
| Too many hits | E-value too lenient | Use `-E 1e-10` or stricter |
| No hits found | HMM profile not calibrated | Try jackhmmer or lower E-value |
| Missing known members | Coverage filter too strict | Check domain coverage |
| Slow search | Large database | Use `--noali` to skip alignment |
| Memory error | Very large database | Split database into chunks |

## Integration with Gene Family Workflow

```bash
# Complete gene family identification workflow

# 1. Get HMM profile
wget https://pfam.xfam.org/family/PF03110/hmm -O SBP.hmm

# 2. HMM search
hmmsearch --cpu 8 -E 1e-5 \
  --tblout hmm_results.tbl \
  --domtblout domain_results.tbl \
  SBP.hmm genome_proteins.fasta

# 3. Extract and validate hits
grep -v "^#" hmm_results.tbl | \
  awk '$5 < 1e-5 {print $1}' | sort -u > family_ids.txt

# 4. Check domain coverage
grep -v "^#" domain_results.tbl | \
  awk '{split($12,a,"-"); split($13,b,"-"); cov=(a[2]-a[1])/$2; if(cov>0.5) print $1}' | \
  sort -u > filtered_ids.txt

# 5. Extract sequences
seqkit grep -f filtered_ids.txt genome_proteins.fasta > family_proteins.fasta
```

## References

- Eddy (2011) Accelerated Profile HMM Searches
- Finn et al. (2016) The Pfam protein families database
- HMMER Documentation: http://hmmer.org/documentation.html

## Related Tools

- **BLAST**: Alternative homology search
- **Pfam**: HMM profile database
- **CD-HIT**: Remove redundant sequences
- **MAFFT**: Multiple sequence alignment for HMM building
