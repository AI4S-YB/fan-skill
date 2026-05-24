# Off-Target Analysis (Cas-OFFinder)

**Goal:** Identify potential off-target sites for candidate sgRNAs across the whole plant genome
**Best for:** Selecting sgRNAs with minimal off-target activity

## Prerequisites

- Cas-OFFinder (http://www.rgenome.net/cas-offinder/)
- Reference genome FASTA
- Candidate sgRNA sequences (20-mer protospacer + NGG PAM)

## Cas-OFFinder Usage

### Web Interface

1. Go to http://www.rgenome.net/cas-offinder/
2. Select species (pre-loaded genomes for major crops)
3. Upload or paste sgRNA sequences (format: 20mer + NRG)
4. Set parameters:
   - Mismatch number: up to 3
   - DNA bulge size: up to 1
   - RNA bulge size: up to 1
5. Submit and download results

### Command Line (Local)

```bash
# Download and compile Cas-OFFinder
git clone https://github.com/snugel/cas-offinder.git
cd cas-offinder
make

# Prepare input file
cat > input.txt << EOF
/path/to/genome.fa
NNNNNNNNNNNNNNNNNNNNNGG 5
sgRNA1 NNNNNNNNNNNNNNNNNNNNNGG
sgRNA2 NNNNNNNNNNNNNNNNNNNNNGG
EOF

# Run Cas-OFFinder
./cas-offinder input.txt C output.txt

# 5 = max mismatches (including DNA/RNA bulge)
```

### Interpreting Results

```
Output format:
sgRNA_name  Chromosome  Position  Target_Sequence  Mismatches  Strand

Example:
sgRNA1  Chr01  12345678  NNNNNNNNNNNNNNNNNGGG  3  +
sgRNA1  Chr02  87654321  NNNNNNNNNANNNNNNNNGG  1  -

Filter:
- 0-1 mismatch off-targets: HIGH RISK — consider rejecting sgRNA
- 2-3 mismatch off-targets: MODERATE RISK — check genomic context
- >3 mismatches: LOW RISK — generally acceptable
```

## Off-Target Filtering Pipeline

```python
import pandas as pd

def filter_off_targets(offtarget_file, max_mismatch=3):
    """Filter and rank sgRNAs by off-target profile"""
    df = pd.read_csv(offtarget_file, sep="\t", header=None,
                     names=["sgRNA", "Chr", "Pos", "Seq", "Mismatch", "Strand"])

    # Count off-targets by mismatch count
    summary = df.groupby(["sgRNA", "Mismatch"])["Seq"].count().unstack(fill_value=0)

    # Scoring: penalty for each off-target
    penalty_weights = {
        0: 100,   # Perfect match (very bad)
        1: 50,    # 1 mismatch (bad)
        2: 10,    # 2 mismatches (moderate)
        3: 1,     # 3 mismatches (minor)
    }

    score = 0
    for mm, weight in penalty_weights.items():
        if mm in summary.columns:
            score += summary[mm] * weight

    # Genic vs intergenic classification (requires GFF)
    # ...

    return summary, score

# Rank sgRNAs by off-target risk
ranked = sorted(off_target_scores.items(), key=lambda x: x[1])
print("Best sgRNAs (lowest off-target risk):")
for sgrna, score in ranked[:5]:
    print(f"  {sgrna}: off-target score = {score}")
```

## CFD Score Calculation

```python
# Cutting Frequency Determination (CFD) score
# Estimates the likelihood of off-target cleavage

def calculate_cfd(mismatch_positions):
    """
    CFD score based on mismatch positions in protospacer
    Reference: Doench et al. 2016, Nature Biotechnology
    """
    # Position-specific weights (1-20 from PAM-distal to PAM-proximal)
    # Higher weight = more tolerance for mismatch at this position
    weights = {
        20: 0.0,  # PAM-proximal (position 20) — least tolerant
        19: 0.0,
        18: 0.1,
        17: 0.2,
        16: 0.3,
        # ...
        1: 0.9,   # PAM-distal (position 1) — most tolerant
    }

    cfd = 1.0
    for pos in mismatch_positions:
        cfd *= (1 - weights.get(pos, 0.5))

    return cfd
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| Max mismatches | 3 | Beyond 3 mismatches, cleavage probability is negligible |
| DNA bulge | 1 | Allow 1-base DNA bulge (RNA-guided indel formation) |
| RNA bulge | 1 | Allow 1-base RNA bulge |
| Genic filter | Exclude CDS off-targets | Off-targets in exons are high-risk |
| CFD threshold | > 0.8 | High specificity |

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| Max mismatches | 3 | For large plant genomes (>1Gb), consider 4 mismatches; for small genomes (<200Mb), limit to 2 | Large genomes have exponentially more chance matches — 3 mismatches in wheat may return thousands of hits; conversely, small genomes need stricter filtering |
| DNA bulge size | 1 | For Cas12a, set to 0 (Cas12a has different indel tolerance); for SpCas9-NG, keep at 1 | Different Cas enzymes have different bulge tolerance profiles; Cas12a is less tolerant of DNA bulges |
| RNA bulge size | 1 | Set to 0 for stringent specificity screening | RNA bulges are biologically possible but rare; removing them reduces false-positive off-target calls |
| Genomic search scope | Nuclear genome only | **Always include chloroplast and mitochondrial genomes** | Cas9/gRNA can enter organelles in plant cells; organellar off-target edits have been experimentally confirmed |
| CFD threshold | >0.8 | For essential/cell-cycle genes, raise to >0.9; for non-coding targets, accept >0.7 | Higher stringency needed when off-target editing in essential genes could be lethal |

## Plant-Specific Notes

- **Large genome scan strategy**: For species with large genomes (wheat 17Gb, barley 5Gb, maize 2.4Gb), run Cas-OFFinder chromosome-by-chromosome rather than on the whole genome. A whole-genome scan with 3 mismatches on wheat can run for days and produce millions of candidate off-targets. Split the genome file by chromosome and run parallel jobs, then merge results.
- Repetitive genome regions: sgRNAs matching repeat elements = hundreds of off-targets -> reject
- Chloroplast/mitochondrial DNA: Cas9 can enter organelles — include organellar genomes in off-target search
- Homeologous genes in polyploids: "off-target" may be intentional (targeting all homeologs)
- Check off-targets in promoter regions (<1kb upstream) — potential regulation disruption
- **Transgene-free editing**: When generating transgene-free edited plants, off-target mutations are heritable and cannot be segregated away. Increase specificity stringency (require specificity score >95, 0 off-targets with <=2 mismatches) for any sgRNA used in transgene-free pipelines.

## Off-Target Risk Classification

```python
def classify_off_target_risk(sgrna_name, off_targets_df):
    """Classify sgRNA off-target risk level"""
    ot = off_targets_df[off_targets_df["sgRNA"] == sgrna_name]

    mm0 = len(ot[ot["Mismatch"] == 0])
    mm1 = len(ot[ot["Mismatch"] == 1])
    mm2 = len(ot[ot["Mismatch"] <= 2])

    if mm1 > 0:
        return "HIGH_RISK"  # Even 1 off-target with 0-1 mismatch = reject
    elif mm2 > 5:
        return "MODERATE_RISK"  # Multiple 2-mm off-targets = caution
    elif mm2 > 0:
        return "LOW_RISK"  # Few 2-mm off-targets = acceptable
    else:
        return "MINIMAL_RISK"  # No off-targets with <=2 mismatches
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No output" from Cas-OFFinder | Input format incorrect | Check: 20mer + 'N' + 'GG' (23-mer total) |
| Thousands of off-targets | sgRNA in repetitive region | Reject sgRNA, design new one |
| Cas-OFFinder crashed | Genome too large | Split genome by chromosome |
| Missing organellar off-targets | Only nuclear genome searched | Append chloroplast/mitochondrial genome |
