# PLINK2 — Genotype Data Management & QC Pipeline

**Goal:** Full genotype data pipeline: format conversion, quality control, LD pruning, PCA, and export for downstream GWAS / population genetics / genomic selection
**Approach:** Convert → QC filter (mind/geno/maf/hwe) → LD prune → PCA → export
**Best for:** Any plant genotype analysis requiring clean, linkage-pruned SNP data

## Prerequisites
- PLINK 1.9 or PLINK 2.0 (`plink` or `plink2` on PATH or in container)
- Input: VCF (gzip or plain), Hapmap, or existing PLINK binary
- For PLINK 1.9: `.ped` + `.map`; for PLINK 2.0: `.pgen` + `.pvar` + `.psam`

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data and analysis type.

```bash
# ── Step 1: Convert VCF to PLINK binary ──
plink2 \
  --vcf ${INPUT}.vcf.gz \
  --allow-extra-chr \
  --make-bed \
  --out ${OUTPUT_PREFIX}_raw

# ── Step 2: Sample QC (filter by missingness) ──
plink2 \
  --bfile ${OUTPUT_PREFIX}_raw \
  --allow-extra-chr \
  --mind ${MIND} \               # 0.1 (GWAS) / 0.2 (pop) / 0.05 (GS)
  --out ${OUTPUT_PREFIX}_sample_qc

# ── Step 3: SNP QC (missingness, MAF, HWE) ──
plink2 \
  --bfile ${OUTPUT_PREFIX}_sample_qc \
  --allow-extra-chr \
  --geno ${GENO} \               # 0.1 (GWAS) / 0.2 (pop) / 0.05 (GS)
  --maf ${MAF} \                 # 0.05 (GWAS) / 0.01 (pop) / 0.05 (GS)
  --hwe ${HWE} \                 # 1e-6 (GWAS) / 1e-50 (pop) / --hwe skip for inbreds
  --out ${OUTPUT_PREFIX}_snp_qc

# ── Step 4: LD pruning ──
plink2 \
  --bfile ${OUTPUT_PREFIX}_snp_qc \
  --allow-extra-chr \
  --indep-pairwise ${WINDOW} ${STEP} ${R2} \   # 50 5 0.2 (standard)
  --out ${OUTPUT_PREFIX}_prune_list

plink2 \
  --bfile ${OUTPUT_PREFIX}_snp_qc \
  --allow-extra-chr \
  --extract ${OUTPUT_PREFIX}_prune_list.prune.in \
  --make-bed \
  --out ${OUTPUT_PREFIX}_pruned

# ── Step 5: PCA ──
plink2 \
  --bfile ${OUTPUT_PREFIX}_pruned \
  --allow-extra-chr \
  --pca ${N_PCS} \
  --out ${OUTPUT_PREFIX}_pca

# ── Step 6: Export formats for downstream tools ──
# VCF (for bcftools, BEAGLE imputation)
plink2 --bfile ${OUTPUT_PREFIX}_pruned --allow-extra-chr --recode vcf --out ${OUTPUT}_vcf

# 012 genotype matrix (for GAPIT, rrBLUP, rrBLUP GS)
plink2 --bfile ${OUTPUT_PREFIX}_pruned --allow-extra-chr --recode A --out ${OUTPUT}_012

# tped/tfam (for EMMAX, TASSEL)
plink2 --bfile ${OUTPUT_PREFIX}_pruned --allow-extra-chr --recode 12 --out ${OUTPUT}_tped
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `--mind` | 0.1 | GWAS → 0.1; Population → 0.2; GS → 0.05 | Population studies tolerate more missing data to retain rare alleles; GS needs clean training data |
| `--geno` | 0.1 | GWAS → 0.1; Population → 0.2; WGS/imputed → 0.05 | Imputed/WGS has more SNPs so stricter filter is affordable |
| `--maf` | 0.05 | GWAS → 0.05; Population → 0.01; Diversity panel → 0.01 | Population structure needs rare variants; GWAS discards them (power too low) |
| `--hwe` | 1e-6 | GWAS → 1e-6; Population → 1e-50; **Inbred/selfing → skip** | HWE is invalid for inbred lines (see Plant-Specific) |
| `--indep-pairwise` | 50 5 0.2 | PCA only → 50 5 0.1; Haplotype blocks → 1000 100 0.8 | Tighter r^2 for PCA independence; looser for retaining LD block signal |
| `--pca` | 10 | Population structure → 10; GWAS covariate → 3-5; GS → 5 | More PCs for fine population structure; fewer for GWAS to avoid over-correction |
| `--double-id` | off | On when VCF has duplicate FID/IID patterns | PLINK sets FID=IID by default; `--double-id` preserves original |
| `--allow-extra-chr` | **always on** | Plants: **always required** | Non-standard chromosome names (Chr01, scaffold_123, NC_xxx) |

### QC Thresholds by Analysis Type

| Filter | GWAS (GEMMA/FarmCPU/BLINK) | Population genetics | Genomic Selection |
|--------|:---:|:---:|:---:|
| `--mind` | 0.1 | 0.2 | 0.05 |
| `--geno` | 0.1 | 0.2 | 0.05 |
| `--maf` | 0.05 | 0.01 | 0.05 |
| `--hwe` | 1e-6 (skip if inbred) | **skip for landraces** | skip for inbred lines |
| LD prune | yes (for PCA) | yes (for PCA/STRUCTURE) | optional |
| `--pca` | 3-5 PCs | 10 PCs | 5 PCs |

---

## Format Conversion Quick Reference

```bash
# VCF → PLINK binary
plink2 --vcf input.vcf.gz --allow-extra-chr --make-bed --out output

# VCF with dosage (DS field) → PLINK
plink2 --vcf input.vcf.gz dosage=DS --allow-extra-chr --make-bed --out output

# PLINK binary → VCF
plink2 --bfile input --allow-extra-chr --recode vcf --out output

# Hapmap → PLINK binary
plink2 --file input --allow-extra-chr --make-bed --out output

# VCF with duplicate variant IDs → regenerate
plink2 --vcf input.vcf.gz --allow-extra-chr --set-all-var-ids @:# --make-bed --out output
```

---

## Plant-Specific Notes

### Chromosome naming chaos (universal)
- Plants use non-standard chromosome names: `Chr01`, `chr1`, `1`, `NC_039463.1`, `scaffold_1234`
- **Always use `--allow-extra-chr`** — without it, PLINK rejects anything not `1-22, X, Y`
- For reference genome scaffolds: expect many "Invalid chromosome code" errors without this flag

### Inbred / selfing species — skip HWE
- Rice, soybean, wheat, barley, tomato, and many other crops are self-pollinated
- **HWE filtering is biologically wrong for inbred lines** — they are homozygous by design
- Use `--hwe 1e-50` as a very lenient filter (catches only genotyping errors), or skip entirely
- Outcrossing species (maize, sunflower, apple, pear): standard HWE filter is appropriate

### Polyploid species (wheat, potato, canola, cotton, sugarcane)
- **PLINK natively handles only diploids** — polyploids require preprocessing
- Tetraploids: code as two separate diploid calls (A/B subgenomes in wheat)
- Hexaploid wheat: split by subgenome (A, B, D) before PLINK
- Alternatively: use `--vcf` + GAPIT (which accepts dosage) instead of PLINK binary format
- Potato (autotetraploid): use `polyRAD` or `updog` for genotype calling, then `--recode A` (dosage)

### Duplicate variant IDs
- Plant VCFs often have duplicate IDs from multi-sample joint calling
- Fix with: `--set-all-var-ids @:#` (chr:pos format) or `--set-all-var-ids @:#_\$r` (chr:pos_ref)
- Check first: `zgrep -v "^#" input.vcf.gz | cut -f3 | sort | uniq -d | head`

### Large genomes / high SNP density
- If `.bim` file is enormous (>2M SNPs), LD prune in chunks: split by chromosome, prune each, then merge
- For WGS data from large genomes (wheat ~17 Gbp, maize ~2.4 Gbp): use `--thin` to downsample first
- PLINK2 handles larger datasets better than PLINK 1.9 — prefer `plink2` for >1M SNPs

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Invalid chromosome code" | Non-standard chr names (Chr01, scaffold_123) | Add `--allow-extra-chr` (plants always need this) |
| "Error: .bim file not found" | Wrong file prefix or files in different directory | Verify with `ls prefix.*`; use full path or cd to file directory |
| "Error: .fam file not found" | Missing `.fam` from incomplete conversion | Re-run `--make-bed`; check input VCF is not truncated |
| "Duplicate variant ID" | VCF has duplicate IDs from multi-sample calling | Add `--set-all-var-ids @:#` to regenerate unique IDs |
| "VCF file not recognized" | VCF version or encoding issue | Check header with `zgrep "^##" file.vcf.gz`; ensure VCFv4.1+ |
| PLINK2 segfaults on large VCF | Memory exhaustion | Split VCF by chromosome; process in parallel; use PLINK2 not PLINK1.9 |
| "Warning: At least one allele is not A/C/G/T" | Indels or structural variants present | Add `--snps-only` to exclude non-SNP variants; or process with `--allow-no-sex` |
| PCA shows batch effect instead of population structure | Unpruned LD or platform batch | LD-prune before PCA; check if genotypes from different sequencing runs |
| QC removes all SNPs | `--maf` or `--geno` too strict | Relax: start with `--geno 0.3 --maf 0.01`, then progressively tighten |
| "Error: .pvar file not found" (PLINK2) | Using PLINK1.9 flags with PLINK2 binary files | Use `plink` for `.bed/.bim/.fam`; use `plink2` with `--pfile` for `.pgen/.pvar/.psam` |
