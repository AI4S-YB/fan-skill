# Selection Pressure Analysis (dN/dS)

**Goal:** Detect positive selection (dN/dS > 1) or purifying selection (dN/dS < 1)
**Best for:** Identifying genes under adaptive evolution

## Prerequisites
- PAML 4.9+ (codeml)
- Codon-aligned CDS sequences
- Phylogenetic tree (Newick format)

## Branch-Site Model

```
PAML control file (codeml.ctl):
model = 2 (branch-site), NSsites = 2
LRT: 2*(lnL_alt - lnL_null) ~ χ²(1)
```

## Plant-Specific Notes
- Plant resistance genes (NBS-LRR) frequently show dN/dS > 1
- Photosynthesis genes show strong purifying selection (dN/dS < 0.1)
- After polyploidization, duplicated genes may show relaxed selection

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "codon alignment error" | Frameshift or stop codon | Check CDS quality |
