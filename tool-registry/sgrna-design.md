# sgRNA Design (CRISPOR / CRISPR-P / CHOPCHOP)

**Goal:** Design optimal sgRNAs for CRISPR/Cas9 genome editing in plant species
**Best for:** Selecting high-efficiency, high-specificity guide RNAs targeting plant genes

## Prerequisites

- Target gene sequence (FASTA or genomic coordinates)
- Species reference genome (available for major crops in CRISPOR)
- Web browser (CRISPOR, CRISPR-P) or Python/Perl (CHOPCHOP local)

## CRISPOR (Major Crop Species)

### Web Interface

1. Go to http://crispor.tefor.net/
2. Select species genome: *Oryza sativa*, *Zea mays*, *Arabidopsis thaliana*, *Glycine max*, *Solanum lycopersicum*, *Triticum aestivum*
3. Enter target sequence (FASTA) or genomic coordinates
4. Select PAM: NGG for SpCas9 (default) or alternative
5. Click "Submit"

### Interpreting Results

CRISPOR provides per-sgRNA scores:
- **Specificity score** (MIT): 0-100, higher = fewer off-targets
- **Efficiency score** (Doench 2016): 0-100, predicted on-target activity
- **Off-targets**: Listed by mismatch count (0,1,2,3,4 mismatches)

### Selection Criteria

```
Priority order:
1. Specificity score > 80
2. Efficiency score > 50
3. 0 off-targets with 0-1 mismatches
4. < 5 off-targets with 0-2 mismatches
5. Target position: early exon (first 50% of CDS)
6. GC content: 40-60%
```

## CRISPR-P (Non-Model Species)

### Web Interface

1. Go to http://crispr.hzau.edu.cn/CRISPR2/
2. Upload custom genome FASTA (if species not listed)
3. Enter gene ID or target sequence
4. Select Cas type: SpCas9, Cas12a, or base editor
5. Select PAM type
6. Submit

### Custom Genome Upload

```bash
# Prepare genome for CRISPR-P
# 1. Concatenate all chromosomes into one FASTA
cat Chr*.fa > genome.fa

# 2. If genome is very large (>3Gb), split by chromosome
# CRISPR-P has file size limits

# 3. Upload via web interface
```

## CHOPCHOP (Local/Command Line)

```bash
# Install
git clone https://bitbucket.org/valenlab/chopchop.git
cd chopchop
pip install -r requirements.txt

# Basic usage
python chopchop.py \
  -G genome.fa \
  -o output_dir/ \
  -t WHOLE \
  -T 1 \
  --scoringMethod DOENCH_2016 \
  -f target_gene.fa

# With genomic coordinates
python chopchop.py \
  -G genome.fa \
  -g genome.gff \
  -o output_dir/ \
  -t CDS \
  -T 1 \
  --scoringMethod DOENCH_2016 \
  --target CHR01:1000000-1005000
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| sgRNA length | 20 nt (SpCas9) | Standard protospacer length |
| PAM | NGG (SpCas9) | Canonical Cas9 PAM |
| GC content | 40-60% | Optimal for sgRNA stability and activity |
| Position | Early CDS (first 30-50%) | Early STOP = complete KO |
| Avoid | TTTT (>=4T) | Pol III terminator signal |
| Avoid | GCG at position 1 | Reduced expression from U6 promoter |

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| Specificity score (MIT) | >80 | For polyploid genomes, require >90; for non-coding targets, accept >70 | Polyploids have homeologous copies that increase off-target risk; non-coding off-targets are less likely to cause functional disruption |
| Efficiency score (Doench) | >50 | For essential genes where partial KO is acceptable, accept >30; for complete KO requirement, require >60 | Lower efficiency sgRNAs may still produce sufficient editing for some applications |
| GC content range | 40-60% | For AT-rich plant genomes (e.g., rice at 43% GC), accept 30-70%; for GC-rich regions, accept up to 75% | Many crop genomes deviate from the 40-60% mammalian optimum; sgRNA expression and stability vary with GC content |
| Target position in CDS | First 30-50% of CDS | Target a constitutive exon shared by all isoforms; avoid the first 5-10% if alternative start codons exist | Early STOP codons ensure complete KO; isoform-specific exons may leave functional splice variants |
| Max off-targets with 0-1mm | 0 | For polyploid homeolog targeting, intentional mismatches to non-target subgenomes are acceptable | Distinguish between undesirable off-targets and intentional multi-copy targeting in polyploids |
| PAM type | NGG (SpCas9) | For AT-rich regions lacking NGG, use Cas9-NG (NG PAM) or Cas12a (TTTV PAM) | AT-rich plant genomes may have sparse NGG sites in target regions |

## Plant-Specific Notes

- For polyploids: check if sgRNA targets one or all subgenome copies. **Polyploid target specificity**: Design sgRNAs that either (a) uniquely target one subgenome by exploiting subgenome-specific SNPs in the protospacer, or (b) target all homeologs simultaneously by selecting conserved regions. Each strategy has different design constraints — specificity vs pan-subgenome editing.
- Use U3/U6 promoters optimized for your species (OsU3 for rice, TaU6 for wheat, AtU6-26 for Arabidopsis)
- For multiplex editing: tRNA-gRNA or Csy4 polycistronic systems
- First exon targeting: ensure the target exon is present in all splice variants
- Check target site in all available cultivars/varieties (SNPs in target site = failed editing)
- **Large plant genomes**: In species with large repetitive genomes (wheat 17Gb, maize 2.4Gb), sgRNAs matching repetitive elements will produce hundreds of off-targets. Always mask repetitive regions before sgRNA design using RepeatMasker or species-specific repeat libraries.

## sgRNA Filtering Script

```python
import re

def score_sgrna(seq, genome_fa):
    """Score candidate sgRNAs"""
    scores = {}

    # GC content
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
    scores['gc'] = gc

    # T-stretch check
    if 'TTTT' in seq:
        scores['t_stretch'] = False
    else:
        scores['t_stretch'] = True

    # PAM check
    if seq.endswith('GG'):
        scores['pam'] = 'NGG'
    else:
        scores['pam'] = 'non-canonical'

    # Position (relative to CDS start, requires annotation)
    # This is gene-specific

    # Self-complementarity
    seq_rc = str(Seq(seq).reverse_complement())
    # ...check for strong secondary structure

    return scores

# Select top sgRNAs
candidates = []  # filled with sgRNA sequences and scores
ranked = sorted(candidates,
                key=lambda x: (x['specificity_score'] * 0.6 + x['efficiency_score'] * 0.4),
                reverse=True)
top3 = ranked[:3]
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| No sgRNA found | No NGG PAM in target region | Use Cas9-NG or Cas12a with alternative PAM |
| All sgRNAs have low specificity | Target in conserved domain or repetitive region | Target less conserved region (UTR, intron) |
| Low predicted efficiency | Unfavorable target sequence | Screen more candidates or shift target region |
| CRISPR-P upload fails | Genome file too large | Split by chromosome, process individually |
