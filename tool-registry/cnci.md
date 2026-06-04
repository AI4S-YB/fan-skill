# CNCI: Coding-Non-Coding Index

## Tool Overview

**Tool ID**: `cnci`
**Category**: lncRNA Identification
**Purpose**: Distinguish coding from non-coding transcripts using sequence features
**Source**: https://github.com/www-bioinfo-org/CNCI

## When to Use

- Long non-coding RNA (lncRNA) identification
- De novo transcriptome annotation
- Distinguishing ncRNA from mRNA
- Filter coding potential in RNA-seq data

## Key Parameters

| Parameter | Default | Recommended for Plants | Description |
|-----------|---------|------------------------|-------------|
| `-m` | ve | **pl** (plants) | Model type |
| `-f` | fasta | fasta | Input file format |
| `-o` | . | output_dir | Output directory |

### Critical Parameter: `-m` (Model Type)

**Plant-specific: Use `-m pl`**

| Model | Description | Applicable To |
|-------|-------------|---------------|
| `ve` | Vertebrate | Animals, human |
| `pl` | **Plant** | Arabidopsis, rice, maize, etc. |
| `mx` | Mixed | Unknown species |

**Why plant model matters**:
- Plant transcripts have different codon usage patterns
- GC content differs from animals
- Plant lncRNAs have distinct sequence features

## Command Examples

### Basic Plant lncRNA Analysis

```bash
# Analyze transcript sequences for coding potential
python CNCI.py -i transcripts.fasta \
  -m pl \
  -f fasta \
  -o cnci_output
```

### Combined with Other Tools

For comprehensive lncRNA identification, combine CPC2 and CNCI:

```bash
# Step 1: Run CPC2
python CPC2.py -i transcripts.fasta -o cpc2_output

# Step 2: Run CNCI
python CNCI.py -i transcripts.fasta -m pl -o cnci_output

# Step 3: Intersect results
# Non-coding = CPC2 non-coding AND CNCI non-coding
```

## Output Interpretation

### CNCI Score

| Score Range | Classification |
|-------------|----------------|
| **> 0** | Coding (positive score) |
| **< 0** | Non-coding (negative score) |
| **≈ 0** | Uncertain |

### Output Format

```
Transcript_ID	CNCI_score	Classification
TRINITY_DN001_c0_g1	-1.234	noncoding
TRINITY_DN002_c0_g1	0.567	coding
TRINITY_DN003_c0_g1	-2.145	noncoding
```

## Plant-Specific Workflow

### Complete lncRNA Pipeline

```bash
# 1. Assemble transcripts (Trinity)
Trinity --seqType fq --left R1.fq --right R2.fq \
  --max_memory 100G --output trinity_out

# 2. Filter by length (>= 200bp)
seqkit seq -m 200 trinity_out/Trinity.fasta > transcripts_filt.fasta

# 3. Predict ORFs (TransDecoder)
TransDecoder.LongOrfs -t transcripts_filt.fasta
TransDecoder.Predict -t transcripts_filt.fasta

# 4. Filter known proteins (BLAST)
blastp -query longest_orfs.pep -db swissprot -out blasp.out

# 5. Remove transcripts with significant protein hits
# Keep transcripts without protein hits

# 6. Run CNCI for remaining transcripts
python CNCI.py -i candidate_lncrnas.fasta \
  -m pl \
  -o cnci_output

# 7. Run CPC2 for comparison
python CPC2.py -i candidate_lncrnas.fasta -o cpc2_output

# 8. Identify high-confidence lncRNAs
# Both CNCI and CPC2 predict non-coding
```

### lncRNA Classification

```bash
# Classify lncRNAs by genomic context
# Requires genome annotation

# lincRNA: intergenic
# intronic lncRNA: within intron
# antisense lncRNA: antisense to gene
```

## Integration with Other Tools

### CPC2 (Coding Potential Calculator 2)

```bash
# CPC2 uses ORF features and isORF
python CPC2.py -i transcripts.fasta -o cpc2_output
```

### FEELnc

```bash
# FEELnc for plant lncRNA
FEELnc_filter.pl -i transcripts.fasta --monoex=-1
```

### Combined Filtering Strategy

```r
# R code to combine results
cnci <- read.table("cnci_output/cnci_result.txt", header=TRUE)
cpc2 <- read.table("cpc2_output/cpc2_result.txt", header=TRUE)

# Merge results
merged <- merge(cnci, cpc2, by="Transcript_ID")

# High-confidence lncRNAs: both tools predict non-coding
high_conf_lncrna <- merged[merged$CNCI_score < 0 & merged$CPC2_score < 0, ]
```

## Plant-Specific Considerations

### 1. Multi-exonic vs Mono-exonic

Plant lncRNAs can be:
- **Multi-exonic**: Spliced transcripts (like mRNA)
- **Mono-exonic**: Single-exon transcripts

Some tools have bias against mono-exonic transcripts. CNCI handles both.

### 2. ORF Length

Plant lncRNAs may have short ORFs (<100aa) but still be non-coding. CNCI considers:
- ORF length
- Codon usage bias
- Sequence composition

### 3. Expression Level

lncRNA expression is typically lower than mRNA. Consider:
- Filter low-expression transcripts first
- Validate lncRNA candidates by expression level

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| Wrong model selected | Default is vertebrate | Use `-m pl` for plants |
| All transcripts classified as coding | Contamination | Check for protein contamination |
| All transcripts non-coding | Short sequences | Ensure sequences are full-length |
| Memory error | Large input file | Split input into batches |

## References

- Sun et al. (2013) CNCI: an effective method to identify noncoding RNAs in plants
- GitHub: https://github.com/www-bioinfo-org/CNCI

## Related Tools

- **CPC2**: Coding Potential Calculator 2
- **FEELnc**: Flexible Extraction and Annotation of LncRNAs
- **CPAT**: Coding Potential Assessment Tool
- **PLEK**: Predictor of Long non-coding RNAs and Messenger RNAs
