# MatrixEQTL — Fast eQTL Mapping with Linear Models

**Goal:** Identify genetic variants associated with gene expression using efficient
linear model testing.

**Best for:** cis-eQTL and trans-eQTL discovery in plant populations with >= 30 samples.

## Prerequisites

- R 4.0+, MatrixEQTL package
- Genotype: PLINK binary or CSV format with SNP x sample matrix
- Expression: gene x sample matrix (normalized, e.g., TPM or RPKM)
- Covariates: sample x covariate matrix (PCA, PEER factors, known confounders)

## Basic Usage

```r
library(MatrixEQTL)

# --- Load data ---
# Genotype: SNP x sample matrix (0/1/2 dosage)
snps <- SlicedData$new()
snps$fileDelimiter <- "\t"
snps$fileOmitCharacters <- "NA"
snps$fileSkipRows <- 1
snps$fileSkipColumns <- 1
snps$fileSliceSize <- 2000
snps$LoadFile("genotype_matrix.txt")

# Expression: gene x sample matrix (log2-transformed)
gene <- SlicedData$new()
gene$fileDelimiter <- "\t"
gene$fileOmitCharacters <- "NA"
gene$fileSkipRows <- 1
gene$fileSkipColumns <- 1
gene$fileSliceSize <- 2000
gene$LoadFile("expression_matrix.txt")

# Covariates: sample x covariate matrix
cvrt <- SlicedData$new()
cvrt$fileDelimiter <- "\t"
cvrt$fileOmitCharacters <- "NA"
cvrt$fileSkipRows <- 1
cvrt$fileSkipColumns <- 1
cvrt$LoadFile("covariates.txt")

# --- Gene/SNP location files ---
# Format: gene_id, chr, start, end
gene_loc <- read.table("gene_location.txt", header = TRUE,
                       stringsAsFactors = FALSE)
# Format: snp_id, chr, position
snp_loc <- read.table("snp_location.txt", header = TRUE,
                      stringsAsFactors = FALSE)

# --- Run eQTL analysis ---
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = "cis_eqtl_results.txt",
  output_file_name.cis = "cis_eqtl_results.txt",
  pvOutputThreshold = 1e-5,        # trans threshold
  pvOutputThreshold.cis = 1e-3,    # cis threshold (relaxed, post-hoc FDR)
  useModel = modelLINEAR,          # linear additive model
  errorCovariance = numeric(),     # numeric() = independent errors
  verbose = TRUE,
  snpspos = snp_loc,
  genepos = gene_loc,
  cisDist = 1e6,                   # 1 Mb cis window
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE
)

# --- Summarize results ---
cat("cis-eQTL tests:", me$cis$ntests, "\n")
cat("cis-eGenes (nominal p < 0.05):",
    length(unique(me$cis$eqtls$gene)), "\n")
cat("trans-eQTL tests:", me$trans$ntests, "\n")
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| `useModel` | `modelLINEAR` | Additive allele dosage; standard for eQTL |
| `cisDist` | 1e6 (1 Mb) | Standard for plants; adjust for species LD decay |
| `pvOutputThreshold.cis` | 1e-3 | Relaxed pre-filter; correct with FDR post-hoc |
| `pvOutputThreshold` | 1e-5 | Trans threshold; stricter to limit output size |
| `min.pv.by.genesnp` | FALSE | Avoid per-gene-SNP minimum p-value; unnecessary |
| `errorCovariance` | `numeric()` | Independent errors; change if related samples |

## Cis vs Trans Classification

After running MatrixEQTL, classify each SNP-gene pair:

```r
# Add distance annotation
cis_results <- me$cis$eqtls
cis_results$snp_chr <- snp_loc$chr[match(cis_results$snps, snp_loc$snpid)]
cis_results$snp_pos <- snp_loc$pos[match(cis_results$snps, snp_loc$snpid)]
cis_results$gene_chr <- gene_loc$chr[match(cis_results$gene, gene_loc$geneid)]
cis_results$gene_tss <- gene_loc$start[match(cis_results$gene, gene_loc$geneid)]

# Trans = different chromosome OR distance > 1 Mb
cis_results$is_trans <- with(cis_results,
  snp_chr != gene_chr | abs(snp_pos - gene_tss) > 1e6
)

# FDR correction per eGene (most conservative)
cis_eqtl <- subset(cis_results, !is_trans)
cis_eqtl$fdr <- p.adjust(cis_eqtl$pvalue, method = "BH")

# Significant cis-eGenes at FDR < 0.05
cis_egenes <- unique(subset(cis_eqtl, fdr < 0.05)$gene)
cat("cis-eGenes (FDR < 0.05):", length(cis_egenes), "\n")
```

## Multiple Testing Correction

The total number of tests in a full eQTL scan is:
```
n_tests = n_snps_cis * n_genes + n_snps_trans * n_genes
```
For a typical plant eQTL study: ~10,000 genes x ~5,000 cis-SNPs = 50M tests.

**Recommended correction strategy** (see also `decision-matrix.yaml` > `multiple_testing`):

```r
# For cis-eQTL: FDR per gene (gene-level significance)
library(qvalue)
cis_eqtl$qvalue <- qvalue(cis_eqtl$pvalue)$qvalue

# For trans-eQTL: stricter threshold
# Bonferroni for genome-wide trans tests
trans_threshold <- 0.05 / me$trans$ntests

# Permutation-based FDR for robust eGene calling
# (separate script — see fastQTL approach)
```

## Plant-Specific Notes

- **Inbred species**: Homozygous genotypes (0/2 coding, no heterozygotes).
  MatrixEQTL additive model works identically; effect size is homozygous difference.
- **Polyploid species**: Run MatrixEQTL per subgenome with subgenome-specific
  SNP and gene location files. Do NOT merge subgenomes.
- **Low sample count**: If n < 30, reduce `pvOutputThreshold.cis` to 1e-2 to limit
  outputs. Expect few or zero significant results — low power is biological, not
  a tool failure.
- **PEER/covariate inclusion**: Always include expression PCs or PEER factors as
  covariates. For n samples, use n/4 PEER factors (or 5-10 PCs for small n).
- **Gene-level FDR**: Use the most significant SNP per gene, then apply BH correction
  across genes. This is the standard eGene definition.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "MatrixEQTL package not found" | Not installed | `install.packages("MatrixEQTL")` |
| "number of columns mismatch" | Sample ID mismatch between genotype and expression | Verify column order matches across all inputs |
| "subscript out of bounds" | Empty SNP/gene location match | Check SNP ID and gene ID match between data and location files |
| Memory error | Too many tests | Reduce `pvOutputThreshold` or increase `fileSliceSize` |
| "NA/NaN in genotype" | Missing genotype values | Impute with mean or mode before loading |
| QQ plot heavy inflation | Missing covariates or population structure | Add PCs and PEER factors to covariate matrix |
| Zero significant results | Low power or overly strict covariate correction | Check n, reduce covariate count, verify expression variance |
