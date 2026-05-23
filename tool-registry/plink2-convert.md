# PLINK2 Format Conversion

**Goal:** Convert genotype data between formats (VCF → BED, Hapmap → BED)
**Best for:** Any format conversion before GWAS

## Prerequisites
- PLINK 1.9 or PLINK 2.0 (`plink` or `plink2` on PATH or in container)

## VCF to PLINK Binary

```bash
plink2 --vcf input.vcf.gz --make-bed --out output_prefix
```

For VCF with dosage data:
```bash
plink2 --vcf input.vcf.gz dosage=DS --make-bed --out output_prefix
```

## PLINK Binary to VCF

```bash
plink2 --bfile input_prefix --recode vcf --out output_prefix
```

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| --vcf | Input VCF (gzip or plain) |
| --bfile | Input PLINK binary prefix |
| --make-bed | Output PLINK binary |
| --recode vcf | Output VCF |
| --double-id | Keep FID and IID as-is |
| --allow-extra-chr | Allow non-standard chromosome names (common in plants) |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Error: .bim file not found" | Wrong prefix or missing file | Verify with `ls prefix.*` |
| "Invalid chromosome code" | Non-standard chr names | Add `--allow-extra-chr` |
| "VCF file not recognized" | VCF version or encoding issue | Check VCF header with `zgrep "^##" file.vcf.gz` |
| "Duplicate variant ID" | VCF has duplicate IDs | Add `--set-all-var-ids @:#` to regenerate IDs |
