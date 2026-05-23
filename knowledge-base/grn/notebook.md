# GRN

GRN -- from expression matrix to regulatory modules

## TF

### MYB

MYB R2R3-MYB

- ABA JA
- MYB 1R-MYB R3-MYB

### WRKY

WRKY WRKYGQK -- W-box (TTGACC/T) SA JA
WRKY 70+ 7

### NAC

NAC (NAM/ATAF/CUC) -- NACRS
- NAC -- NAC
- NAC

### bHLH

bHLH basic helix-loop-helix JA MYC2

### AP2/ERF

AP2/ERF -- GCC-box DRE/CRT
- DREB CBF DREB1A/CBF3
- ERF ERF1 PDF1.2 JA/ET

### bZIP

bZIP -- ABRE (ABA-responsive element) ABA ABI5 AREB/ABF

### GRAS

GRAS GA DELLA protein -- GA

### MADS-box

MADS-box -- CArG-box floral organ identity ABCDE model

## TC

GENIE3 tree-based

SCENIC GENIE3 +

- **RcisTarget** cis-regulatory motif
- **AUCell** regulon -- AUC
- **RSS** regulon specificity score

WGCNA TF -- WGCNA module hub

## ploidy

- **allopolyploid** homeolog -- A/B/D genome TF
  - homeolog

- **autopolyploid** allele allelic dosage

## non-model

1. homolog inference Arabidopsis thaliana
2. PlantTFDB
3. Top TF hub known TF family

## GO

GO BP -- TF hub DEG

## hub

hub gene top 5% degree centrality + betweenness centrality Z-score > 2
top hub TF check expression correlation with target genes -- regulon -- in silico

## pitfalls

1. **TF annotation** -- co-expression false positive
2. **low sample** -- GENIE3 n >= 15 < 30 --prior
3. **batch effect** -- batch
4. **polyploid mapping** -- reads homeolog multi-mapping -- expression bias quantification
