# Fst (Fixation Index)

**Goal:** Quantify genetic differentiation between populations
**Best for:** Comparing predefined groups (breed types, geographic origins, subpopulations)

## Prerequisites
- PLINK 1.9+
- Population assignment file

## Basic Usage

```bash
plink --bfile input_ld --fst --within pops.txt --out fst_result
```

pops.txt format:
```
FID IID POP
sample1 sample1 popA
sample2 sample2 popB
```

## Plant-Specific Notes
- Weir & Cockerham Fst is preferred for small/unequal sample sizes
- Inbreeding inflates Fst — report along with heterozygosity
- Plant breeding populations: Fst 0.05-0.15 is typical; wild populations: 0.2+

## Interpretation
| Fst Range | Differentiation Level |
|-----------|----------------------|
| < 0.05 | Low |
| 0.05-0.15 | Moderate |
| 0.15-0.25 | Great |
| > 0.25 | Very great |
