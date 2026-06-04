# PlantCARE: Plant Cis-Acting Regulatory Element Analysis

## Tool Overview

**Tool ID**: `plantcare`
**Category**: Regulatory Element Analysis
**Purpose**: Identify cis-acting regulatory elements in plant promoter sequences
**Web Server**: http://bioinformatics.psb.ugent.be/webtools/plantcare/html/

## When to Use

- Gene family promoter analysis
- Transcription factor binding site prediction
- Stress-responsive element identification
- Comparative promoter analysis across species

## Key Parameters

| Parameter | Default | Recommended for Plants | Description |
|-----------|---------|------------------------|-------------|
| Promoter length | 1000bp | **2000bp** | Upstream sequence length from ATG |
| Min score | 5 | 5 | Minimum similarity score for motif match |
| Strand | both | both | Search on both DNA strands |

### Critical Parameter: Promoter Length

**Why 2000bp for plants?**

Plant core promoters are typically more compact than animal promoters, but regulatory elements can be found further upstream:

- **TATA box**: -30 to -100bp (core promoter)
- **CAAT box**: -80 to -300bp (proximal promoter)
- **Stress-responsive elements (ABRE, MBS, etc.)**: -200 to -2000bp
- **Enhancer-like elements**: Can be >2000bp upstream

For comprehensive plant promoter analysis, extract **2000bp upstream** of the translation start site (ATG).

## Command Line Usage

### Step 1: Extract Promoter Sequences

```bash
# Extract 2000bp upstream of each gene
# Using BEDTools

# Create gene regions BED file
awk '$3=="gene" {print $1"\t"$4-1"\t"$5"\t"$9}' annotation.gff3 | \
  sed 's/ID=//;s/;.*//' > genes.bed

# Create promoter regions (2000bp upstream)
# For genes on + strand
awk '$6=="+"' genes.bed | \
  awk '{start=$2-2000; if(start<0) start=0; print $1"\t"start"\t"$2"\t"$4}' > promoters_plus.bed

# For genes on - strand (promoter is downstream of gene end)
awk '$6=="-"' genes.bed | \
  awk '{print $1"\t"$3"\t"$3+2000"\t"$4}' > promoters_minus.bed

# Combine and extract sequences
cat promoters_plus.bed promoters_minus.bed | sort -k1,1 -k2,2n > promoters.bed
bedtools getfasta -fi genome.fasta -bed promoters.bed -fo promoters.fasta
```

### Step 2: Run PlantCARE

PlantCARE is primarily a web-based tool. For batch processing:

```bash
# Option 1: Use local PlantCARE installation (if available)
perl plantcare.pl promoters.fasta output_dir/

# Option 2: Submit to web server manually
# Upload promoters.fasta to http://bioinformatics.psb.ugent.be/webtools/plantcare/html/
```

## Output Interpretation

### Common Plant Cis-Elements

| Element | Function | Sequence Motif | Notes |
|---------|----------|----------------|-------|
| **ABRE** | ABA-responsive | ACGTG(G/T)(C/A) | Abscisic acid response |
| **MBS** | Drought-responsive | CAACTG | MYB binding site |
| **G-box** | Light-responsive | CACGTG | bZIP/bHLH binding |
| **W-box** | Wound/pathogen | TTGAC(C/T) | WRKY binding |
| **HSE** | Heat shock | A(A/T)A(A/T)A(T/C)A(A/T)T | Heat stress |
| **TCA-element** | Salicylic acid | CC(T/A)TCTTTTT | SA response |
| **TGA-element** | Auxin-responsive | AACGAC | Auxin response |
| **CGTCA-motif** | MeJA-responsive | CGTCA | Methyl jasmonate |
| **TGACG-motif** | MeJA-responsive | TGACG | Methyl jasmonate |
| **Skn-1-like** | Endosperm expression | GTCAT | Seed-specific |
| **GCN4_motif** | Endosperm expression | TGTGTAG | Seed-specific |
| **ARE** | Anaerobic induction | AAACCA | Low oxygen |

### Output File Format

```
# SequenceName   ElementName   Start   End   Strand   Sequence   Score
Gene001_promoter  ABRE         156     161   +        ACGTGG     6.2
Gene001_promoter  MBS          423     428   -        CAACTG     5.8
```

## Plant-Specific Considerations

### 1. Multi-gene Families

For gene families (e.g., SPL, NAC, WRKY):
- Extract promoters for all family members
- Compare element distribution across subgroups
- Identify subgroup-specific regulatory patterns

### 2. Polyploid Species

In polyploids (wheat, cotton, canola):
- Homeologs may have diverged promoters
- Compare A/B/D (wheat) or At/Dt (cotton) subgenome elements
- Look for regulatory subfunctionalization

### 3. Tissue-Specific Elements

Plant promoters often contain tissue-specific elements:
- Endosperm: Skn-1-like, GCN4_motif
- Root: ROOTMOTIFTAPOX1
- Leaf/Guard cell: G-box, GT1-motif
- Pollen: POLLEN1LELAT52

## Integration with Other Analyses

### Combine with Expression Data

```r
# Correlate promoter elements with expression patterns
# R code example

# Count ABRE elements per promoter
abre_count <- sapply(promoter_elements, function(x) {
  sum(grepl("ABRE", x$ElementName))
})

# Test correlation with drought-induced expression
cor.test(abre_count, log2FC_drought)
```

### Combine with ChIP-seq

If you have ChIP-seq data for a TF:
- Check if predicted binding sites overlap with ChIP peaks
- Validate predicted elements experimentally

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| No elements found | Promoter extraction failed | Check strand orientation; verify ATG position |
| Too many elements | Low stringency | Increase min score threshold |
| Missing known elements | Promoter too short | Extend to 2000bp or more |
| False positives | Repetitive sequences | Mask repeats before analysis |

## References

- Lescot et al. (2002) PlantCARE, a database of plant cis-acting regulatory elements
- PlantCARE Database: http://bioinformatics.psb.ugent.be/webtools/plantcare/html/

## Related Tools

- **PLACE**: Alternative plant cis-element database
- **AGRIS**: Arabidopsis Gene Regulatory Information Server
- **JASPAR Plants**: Plant TF binding profile database
