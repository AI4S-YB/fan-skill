# Gene Family Identification

**Goal:** Identify and classify gene families — groups of homologous genes
**Best for:** Understanding gene expansion/contraction, functional diversification

## Prerequisites
- Protein sequences (one or multiple species)
- HMMER for domain-based search
- BLASTP + MCL for clustering

## HMM Search

```bash
hmmscan --domtblout output.domtbl Pfam-A.hmm proteins.faa
```

## BLAST + MCL Clustering

```bash
blastp -query proteins.faa -db proteins.faa -out blast.txt -outfmt 6
mcxload -abc blast.txt --stream-mirror --stream-neg-log10
mcl input.mci -I 2.0 -use-tab output.tab
```

## Plant-Specific Notes
- NBS-LRR (disease resistance) is the largest and most dynamic plant gene family
- Use specific tools for plant resistance genes: NLR-Annotator, RGAugury, DRAGO2
- Gene family expansion in plants often correlates with stress adaptation

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "too many clusters" | MCL inflation too low | Increase -I to 2.0-4.0 |
