# eQTL Data Preprocessing

**Goal:** Prepare genotype and expression data for eQTL analysis — QC, normalization,
covariate estimation, and data formatting.

**Best for:** All eQTL studies; critical quality control step before MatrixEQTL or
similar tools.

## Prerequisites

- R 4.0+, data.table, dplyr, peer (optional, for PEER factors)
- Genotype data in VCF or PLINK format
- Expression matrix (genes x samples, normalized counts or TPM)
- Sample metadata (IDs linking genotype and expression)

## Workflow

### Step 1: Genotype QC and Filtering

```r
# After converting VCF to PLINK binary (see GWAS Skill):
# plink2 --vcf genotypes.vcf.gz --make-bed --out eqtl_geno

# In R: load PLINK-format genotype data
library(data.table)

# Load .raw file (additive dosage: 0/1/2)
geno <- fread("eqtl_geno.raw", data.table = FALSE)
rownames(geno) <- geno$IID
geno <- geno[, -(1:6)]  # Remove PLINK metadata columns

# Filter by MAF
maf <- colMeans(geno, na.rm = TRUE) / 2
maf <- pmin(maf, 1 - maf)  # Minor allele frequency
geno <- geno[, maf >= 0.05]

# Filter by call rate (SNP level)
snp_call_rate <- 1 - colMeans(is.na(geno))
geno <- geno[, snp_call_rate >= 0.9]

# Impute missing genotypes with mean
for (j in seq_len(ncol(geno))) {
  geno[is.na(geno[, j]), j] <- mean(geno[, j], na.rm = TRUE)
}

cat("SNPs after QC:", ncol(geno), "\n")
cat("Samples:", nrow(geno), "\n")
```

### Step 2: Expression QC and Normalization

```r
# Load expression matrix (genes x samples)
expr <- fread("expression_matrix.txt", data.table = FALSE)
rownames(expr) <- expr$gene_id
expr <- expr[, -1]

# Filter lowly expressed genes (>= 1 TPM in >= 20% of samples)
expr_tpm <- expr  # Assuming data is already in TPM
gene_expr_rate <- rowMeans(expr_tpm >= 1)
expr_filtered <- expr_tpm[gene_expr_rate >= 0.2, ]

cat("Genes after filtering:", nrow(expr_filtered), "\n")

# Log2 transformation
expr_log2 <- log2(expr_filtered + 1)

# Quantile normalization (optional, between-sample)
library(preprocessCore)
expr_normalized <- normalize.quantiles(as.matrix(expr_log2))
rownames(expr_normalized) <- rownames(expr_log2)
colnames(expr_normalized) <- colnames(expr_log2)

cat("Expression matrix:", nrow(expr_normalized), "genes x",
    ncol(expr_normalized), "samples\n")
```

### Step 3: Match Samples Between Genotype and Expression

```r
# Find shared samples
geno_samples <- rownames(geno)
expr_samples <- colnames(expr_normalized)
shared_samples <- intersect(geno_samples, expr_samples)

cat("Genotype samples:", length(geno_samples), "\n")
cat("Expression samples:", length(expr_samples), "\n")
cat("Shared samples:", length(shared_samples), "\n")

if (length(shared_samples) < length(expr_samples) * 0.5) {
  warning("Less than 50% sample overlap — check sample ID format")
}

# Subset to shared samples
geno <- geno[shared_samples, ]
expr_normalized <- expr_normalized[, shared_samples]
```

### Step 4: Covariate Estimation — PEER Factors

```r
# PEER: Probabilistic Estimation of Expression Residuals
# Install: devtools::install_github("PMBio/peer")

# Determine number of PEER factors (rule: n_samples / 4)
n_peer <- floor(ncol(expr_normalized) / 4)
n_peer <- min(n_peer, 30)  # Cap at 30 factors
n_peer <- max(n_peer, 3)   # Minimum 3 factors

cat("Computing", n_peer, "PEER factors...\n")

library(peer)

model <- PEER()
PEER_setPhenoMean(model, as.matrix(t(expr_normalized)))
PEER_setNk(model, n_peer)
PEER_setNmax_iterations(model, 1000)
PEER_update(model)

# Extract factors
peer_factors <- PEER_getX(model)
colnames(peer_factors) <- paste0("PEER", 1:n_peer)
rownames(peer_factors) <- colnames(expr_normalized)

cat("PEER factors shape:", dim(peer_factors), "\n")
```

### Step 5: Covariate Estimation — PCA (Alternative)

```r
# PCA on expression matrix (simpler but less effective than PEER)
expr_centered <- t(scale(t(expr_normalized), scale = FALSE))
pca <- prcomp(t(expr_centered), center = TRUE, scale. = TRUE)

# Number of PCs: variance explained threshold
var_explained <- summary(pca)$importance["Cumulative Proportion", ]
n_pcs <- min(which(var_explained > 0.5), 20)  # Cap at 20 PCs

cat("Using", n_pcs, "PCs (explaining",
    round(var_explained[n_pcs] * 100, 1), "% variance)\n")

pc_covariates <- pca$x[, 1:n_pcs, drop = FALSE]
```

### Step 6: Format Output for MatrixEQTL

```r
# Genotype matrix: transpose to SNP x sample
geno_t <- as.data.frame(t(geno))
geno_t$snp_id <- rownames(geno_t)
geno_t <- geno_t[, c("snp_id", shared_samples)]

write.table(geno_t, "genotype_matrix_eqtl.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Expression matrix: gene x sample
expr_out <- as.data.frame(expr_normalized)
expr_out$gene_id <- rownames(expr_out)
expr_out <- expr_out[, c("gene_id", shared_samples)]

write.table(expr_out, "expression_matrix_eqtl.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Covariate matrix: sample x covariate
covariates <- as.data.frame(t(peer_factors))
covariates$covariate_id <- rownames(covariates)
covariates <- covariates[, c("covariate_id", shared_samples)]

write.table(covariates, "covariates_eqtl.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Files written for MatrixEQTL input.\n")
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| MAF filter | >= 0.05 | Standard; reduce to 0.01 for small populations |
| SNP call rate | >= 0.9 | Remove poorly genotyped variants |
| Gene expression filter | >= 1 TPM in >= 20% samples | Remove unexpressed genes |
| Expression transform | log2(TPM + 1) | Variance stabilization |
| PEER factors: n | n_samples / 4 | Standard in GTEx/eQTL studies |
| PEER factors: max | 30 | Diminishing returns beyond 30 |
| PCA PCs: n | Cumulative variance > 50% | Simpler alternative to PEER |
| Quantile normalization | Optional | Helps with between-sample technical variation |

## Plant-Specific Notes

- **Field-collected samples**: Expression data from field trials has higher technical
  variance. Use more PEER factors (n/3 instead of n/4) and stricter gene expression
  filtering (>= 1 TPM in >= 30% of samples).
- **RNA-seq from different batches**: If samples were sequenced in multiple batches,
  include batch as a known covariate in addition to PEER factors. Batch effects in
  plant studies are common and can completely confound eQTL signals.
- **Tissue-specific filtering**: If analyzing multiple tissues, filter genes
  separately per tissue. A gene silent in leaves may be highly expressed in roots.
  Per-tissue filtering prevents losing tissue-specific eGenes.
- **Ploidy-aware genotype coding**: For polyploids, if using dosage (0/1/2/3/4...),
  MatrixEQTL additive model still works. If using presence/absence, convert to 0/1.
- **Pooled samples**: Some plant studies pool multiple plants per RNA-seq sample.
  This reduces expression variance but also reduces the effective sample size for
  eQTL mapping. Treat each pool as one biological replicate.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Sample IDs don't match" | Different naming conventions between data types | Create a sample mapping file |
| "peer package not found" | PEER not installed | Fall back to PCA (`prcomp` approach) |
| PEER: "system is computationally singular" | Too few samples for requested factors | Reduce n_peer; minimum is 3 |
| Expression variance near zero after filtering | Overly strict gene filter | Relax: >= 0.5 TPM in >= 10% samples |
| Genotype matrix too sparse | Low-coverage sequencing | Impute or use alternative genotyping |
| Covariates remove all eQTL signal | Over-correction (too many PEER factors) | Reduce n_peer; compare with and without |
