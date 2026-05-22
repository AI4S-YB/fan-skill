# Primer3 — PCR Primer Design

**Goal:** Design PCR primers for SNP/InDel/SSR markers
**Best for:** All marker types requiring flanking primers

## Prerequisites
- Primer3-py (Python) or primer3_core (CLI)
- Target sequence in FASTA format

## Basic Usage

```python
import primer3

primers = primer3.design_primers(
    seq_args={
        'SEQUENCE_TEMPLATE': target_seq,
        'SEQUENCE_TARGET': [snp_pos, 1],  # SNP at this position
    },
    global_args={
        'PRIMER_PRODUCT_SIZE_RANGE': [[80, 150]],
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
    }
)
```

## Plant-Specific Notes
- For polyploids: design primers in subgenome-specific regions (SNPs unique to one subgenome)
- For functional markers (dCAPS): introduce mismatches to create restriction sites at the SNP
- KASP markers need two allele-specific forward primers + one common reverse

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "No primers found" | Repeat-rich region | Extend flanking sequence or relax constraints |
