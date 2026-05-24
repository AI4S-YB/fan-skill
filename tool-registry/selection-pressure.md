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

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| model (codeml) | 2 (branch-site) | Use model = 0 for site models testing pervasive selection across all branches; model = 1 for branch models testing selection on specific lineages | Branch-site models (model=2) test for positive selection on specific foreground branches at a subset of sites — most relevant for detecting adaptive evolution in specific plant lineages |
| fix_omega | 0 (alternative); 1 (null) | Always run both: fix_omega=0 (alternative, allows omega>1) and fix_omega=1 (null, omega fixed at 1) | The LRT test requires comparing the alternative model (allows positive selection) against the null model (no positive selection); 2*(lnL_alt - lnL_null) follows chi-square(1) |
| NSsites | 2 | Use NSsites = 0,1,2,7,8 depending on the hypothesis; model 7 (beta) vs 8 (beta+omega) for site-model LRT | Different models test different selection regimes — model 7/8 pair is more conservative than M1a/M2a for site models |
| kappa (ts/tv) | estimated | Fix kappa = 2.0 if estimation fails due to short sequences or low divergence | The transition/transversion ratio must be reasonably estimated for reliable dN/dS inference |
| codonFreq | 2 (F3x4) | Use codonFreq = 3 (F61) for highly biased plant genomes (e.g., grasses with high GC3 content) | Codon frequency model affects dS estimation; F61 uses 61 codon frequencies vs F3x4 uses nucleotide frequencies at 3 codon positions |
| cleandata | 1 | Set cleandata = 0 if removing ambiguous sites deletes most of the alignment | Plant alignments often contain gaps; complete deletion (cleandata=1) may discard informative sites |

## Plant-Specific Notes
- Plant resistance genes (NBS-LRR) frequently show dN/dS > 1, particularly in the LRR domain. **NBS-LRR positive selection**: The solvent-exposed residues of the LRR domain are under diversifying selection in plant R-genes. When analyzing NBS-LRR families, run branch-site models with the LRR domain as foreground to detect specific positively selected sites involved in pathogen recognition.
- Photosynthesis genes (Rubisco, LHC) show strong purifying selection (dN/dS < 0.1)
- After polyploidization, duplicated genes may show relaxed selection (dN/dS approaching 1) before subfunctionalization or neofunctionalization. Compare dN/dS of duplicate pairs pre- and post-duplication to test for selection relaxation.
- Plant transcription factors in stress-responsive families (WRKY, NAC, MYB) often show elevated dN/dS in their activation domains but purifying selection in DNA-binding domains — analyze domains separately.
- Self-incompatibility loci (S-locus) in plants show extreme positive selection (dN/dS often > 2-3) and are positive controls for branch-site model validation.

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "codon alignment error" | Frameshift or stop codon | Check CDS quality; use PAL2NAL for protein-guided codon alignment |
| LRT not significant despite high dN/dS | Too few sites or too few sequences | Increase taxon sampling; run site-model tests (M7 vs M8) which have more power with fewer branches |
| "kappa estimate out of range" | Sequences too divergent or too similar | Fix kappa to 2.0; check if sequences are orthologous |
| "omega > 1 but not significant" | dN/dS inflated by a few sites | Examine BEB posterior probabilities for individual positively selected sites |
| Convergence failure in codeml | Complex model (NSsites=2) with insufficient data | Switch to simpler model (NSsites=0); check for missing data in alignment |
