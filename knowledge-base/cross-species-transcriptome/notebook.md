# Cross-Species Transcriptome Analysis Notebook

## Overview

This notebook provides step-by-step guidance for cross-species transcriptome comparison analysis, covering ortholog identification, expression normalization, batch correction, tissue-specificity analysis, and expression conservation assessment.

---

## Phase 0: Data Preparation

### 0.1 Reference Data Acquisition

```bash
# Download Ensembl reference genomes and annotations (example: human, mouse)
# Using Ensembl REST API or wget
wget ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Download GTF annotations
wget ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
```

### 0.2 Extract Longest Transcript per Gene

```bash
# Use gtftk or custom script to extract longest transcript
# For each species:
gtftk longest_isoform -i annotation.gtf -o longest_transcript.gtf

# Extract CDS and protein sequences
gffread longest_transcript.gtf -g genome.fa -x cds.fa -y protein.fa
```

### 0.3 Expression Data Preparation

```r
# Load expression matrices
human_exp <- read.table("human_tpm.tsv", header=TRUE, row.names=1)
mouse_exp <- read.table("mouse_tpm.tsv", header=TRUE, row.names=1)

# Check data dimensions and distributions
dim(human_exp)  # genes × samples
summary(colSums(human_exp))

# Verify sample metadata
human_meta <- read.table("human_metadata.tsv", header=TRUE)
table(human_meta$tissue)
```

---

## Phase 1: Ortholog Identification

### 1.1 Method Selection

| Method | When to Use | Pros | Cons |
|--------|-------------|------|------|
| OrthoFinder | De novo, multiple species | Accurate, builds species tree | Computationally intensive |
| biomaRt | Ensembl-annotated species | Fast, curated orthologs | Limited to Ensembl species |
| RBH (BLAST) | Quick 1:1 mapping needed | Simple, fast | May miss many-to-many orthologs |
| OrthoDB | Pre-computed orthologs needed | Comprehensive database | May not have latest annotations |

### 1.2 OrthoFinder Pipeline

```bash
# Prepare protein sequence files (one per species)
# Naming convention: SpeciesName.fa

# Run OrthoFinder
orthofinder -f protein_sequences/ \
    -S diamond \
    -M msa \
    -t 32 \
    -a 16

# Output includes:
# - Orthogroups.csv: Gene families
# - Orthogroups/Orthogroups_SingleCopyOrthologues.txt: 1:1 orthologs
# - Species_Tree/SpeciesTree_rooted.txt: Inferred species tree
```

### 1.3 biomaRt for Ensembl Species

```r
library(biomaRt)

# Connect to Ensembl
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get 1:1 orthologs
orthologs <- getLDS(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "with_ortholog",
    values = TRUE,
    mart = human,
    attributesL = c("ensembl_gene_id", "external_gene_name"),
    martL = mouse
)

# Filter for 1:1 orthologs
one_to_one <- orthologs[!duplicated(orthologs$Gene.stable.ID) & 
                        !duplicated(orthologs$Gene.stable.ID.1), ]

# Save mapping
write.table(one_to_one, "human_mouse_orthologs_1to1.tsv", 
            sep="\t", quote=FALSE, row.names=FALSE)
```

### 1.4 Reciprocal Best Hit (RBH) Method

```bash
# Make BLAST databases
makeblastdb -in human_protein.fa -dbtype prot -parse_seqids
makeblastdb -in mouse_protein.fa -dbtype prot -parse_seqids

# Run reciprocal BLAST
blastp -query human_protein.fa -db mouse_protein.fa \
    -evalue 1e-10 -outfmt 6 -num_threads 16 \
    -out human_to_mouse.tsv

blastp -query mouse_protein.fa -db human_protein.fa \
    -evalue 1e-10 -outfmt 6 -num_threads 16 \
    -out mouse_to_human.tsv

# Find RBH pairs (custom script)
# rbh_finder.py human_to_mouse.tsv mouse_to_human.tsv > orthologs_rbh.tsv
```

---

## Phase 2: Expression Normalization

### 2.1 TPM Normalization

```r
# Convert counts to TPM if needed
counts_to_tpm <- function(counts, gene_lengths) {
    # Normalize for gene length (RPK)
    rpk <- counts / (gene_lengths / 1000)
    # Normalize for sequencing depth
    scaling_factor <- colSums(rpk) / 1e6
    tpm <- sweep(rpk, 2, scaling_factor, "/")
    return(tpm)
}

# Log transform with pseudocount
log_tpm <- log2(tpm + 1)
```

### 2.2 Quantile Normalization

```r
library(preprocessCore)

# Apply quantile normalization across species
combined_exp <- cbind(human_log_tpm, mouse_log_tpm)
quantile_normalized <- normalize.quantiles(combined_exp)

# Restore column names
colnames(quantile_normalized) <- colnames(combined_exp)
```

### 2.3 TMM Normalization (for count data)

```r
library(edgeR)

# Create DGEList object
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge, method = "TMM")

# Convert to CPM
cpm <- cpm(dge, log = TRUE, prior.count = 1)
```

---

## Phase 3: Batch Effect Correction

### 3.1 ComBat for Species Batch Correction

```r
library(sva)

# Prepare expression matrix (genes × samples)
# Must use ortholog genes only
ortho_exp <- combined_exp[ortholog_gene_ids, ]

# Create batch variable (species)
batch <- c(rep("human", n_human_samples), 
           rep("mouse", n_mouse_samples))

# Create model matrix for biological variables to preserve
# Example: preserve tissue effects
tissue <- c(human_metadata$tissue, mouse_metadata$tissue)
mod <- model.matrix(~ tissue)

# Apply ComBat
corrected_exp <- ComBat(
    dat = ortho_exp,
    batch = batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
)
```

### 3.2 ComBat-seq (for count data)

```r
library(sva)

# ComBat-seq preserves integer counts
corrected_counts <- ComBat_seq(
    counts = count_matrix,
    batch = batch,
    group = tissue
)
```

### 3.3 MNN Correction (for non-linear effects)

```r
library(batchelor)

# MNN correction for high-dimensional data
corrected_exp <- mnnCorrect(
    human_exp, mouse_exp,
    batch = batch,
    k = 20,  # number of mutual nearest neighbors
    sigma = 0.1  # bandwidth for Gaussian smoothing
)
```

---

## Phase 4: Dimensionality Reduction

### 4.1 PCA Analysis

```r
# Perform PCA
pca_result <- prcomp(
    t(corrected_exp),  # samples in rows
    scale. = TRUE,
    center = TRUE
)

# Variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100

# Visualize
library(ggplot2)
pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    species = batch,
    tissue = tissue
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = species, shape = tissue)) +
    geom_point(size = 3) +
    labs(x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    theme_bw()
```

### 4.2 t-SNE Visualization

```r
library(Rtsne)

# Calculate perplexity based on sample size
n_samples <- ncol(corrected_exp)
perplexity_value <- max(5, floor(sqrt(n_samples)))

# Run t-SNE
set.seed(42)  # for reproducibility
tsne_result <- Rtsne(
    t(corrected_exp),
    dims = 2,
    perplexity = perplexity_value,
    max_iter = 1000,
    check_duplicates = FALSE
)

# Visualize
tsne_df <- data.frame(
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2],
    species = batch,
    tissue = tissue
)

ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = species, shape = tissue)) +
    geom_point(size = 3) +
    theme_bw()
```

### 4.3 UMAP Visualization

```r
library(umap)

# Set UMAP parameters
umap_config <- umap.defaults
umap_config$n_neighbors <- 15
umap_config$min_dist <- 0.1
umap_config$metric <- "euclidean"

# Run UMAP
umap_result <- umap(t(corrected_exp), config = umap_config)

# Visualize
umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    species = batch,
    tissue = tissue
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = species, shape = tissue)) +
    geom_point(size = 3) +
    theme_bw()
```

---

## Phase 5: Tissue Complexity Exploration

### 5.1 Expressed Genes per Tissue

```r
# Define expression threshold
tpm_threshold <- 1  # TPM > 1 considered expressed

# Count expressed genes per tissue/species
expressed_counts <- sapply(unique(tissue), function(t) {
    tissue_samples <- which(tissue == t)
    rowSums(corrected_exp[, tissue_samples] > tpm_threshold) / length(tissue_samples)
})

# Calculate mean expressed genes per tissue
mean_expressed <- apply(expressed_counts, 2, mean)

# Barplot
barplot(mean_expressed, 
        main = "Expressed Genes per Tissue",
        ylab = "Number of Expressed Genes",
        col = "steelblue",
        las = 2)
```

### 5.2 Expression Abundance Distribution

```r
# Distribution of expression values by species
par(mfrow = c(1, 2))

# Human
hist(log_tpm_human, breaks = 100, 
     main = "Human Expression Distribution",
     xlab = "log2(TPM+1)", col = "lightblue")

# Mouse
hist(log_tpm_mouse, breaks = 100,
     main = "Mouse Expression Distribution",
     xlab = "log2(TPM+1)", col = "lightcoral")
```

### 5.3 Dispersion Analysis

```r
# Calculate coefficient of variation per gene
calculate_cv <- function(x) {
    sd(x) / mean(x) * 100
}

cv_human <- apply(corrected_exp[, human_samples], 1, calculate_cv)
cv_mouse <- apply(corrected_exp[, mouse_samples], 1, calculate_cv)

# Identify high-variance genes
high_var_threshold <- quantile(c(cv_human, cv_mouse), 0.9)
high_var_genes <- names(which(cv_human > high_var_threshold | 
                               cv_mouse > high_var_threshold))
```

---

## Phase 6: Tissue-Specificity Analysis

### 6.1 Tau Index Calculation

```r
# Tau tissue-specificity index
# tau = (1 - max(tissue_exp) / sum(tissue_exp)) * (n_tissues - 1) / n_tissues
# Alternative formula: tau = sum(1 - x/max(x)) / (n-1)

calculate_tau <- function(expression_vector) {
    if (max(expression_vector) == 0) return(NA)
    normalized <- expression_vector / max(expression_vector)
    tau <- sum(1 - normalized) / (length(normalized) - 1)
    return(tau)
}

# Calculate tau for each gene in each species
# First, calculate mean expression per tissue
tissue_mean_human <- tapply(corrected_exp[, human_samples], 
                            human_metadata$tissue, 
                            colMeans)

tau_human <- apply(tissue_mean_human, 1, calculate_tau)
tau_mouse <- apply(tissue_mean_mouse, 1, calculate_tau)

# Compare tau between orthologs
tau_comparison <- data.frame(
    gene = ortholog_gene_ids,
    tau_human = tau_human[ortholog_gene_ids],
    tau_mouse = tau_mouse[ortholog_gene_ids]
)

# Correlation of tau values
cor(tau_comparison$tau_human, tau_comparison$tau_mouse, use = "complete.obs")
```

### 6.2 Using tspex Package

```r
library(tspex)

# Prepare tissue mean expression matrix
tissue_means <- aggregate(corrected_exp, 
                         by = list(tissue = tissue), 
                         FUN = mean)
rownames(tissue_means) <- tissue_means$tissue
tissue_means <- t(tissue_means[, -1])

# Calculate tau using tspex
tau_values <- tissue_specificity(tissue_means, method = "tau")

# Calculate SPM (Scaled Specificity Measure)
spm_values <- tissue_specificity(tissue_means, method = "spm")

# Visualize distribution
par(mfrow = c(1, 2))
hist(tau_values, breaks = 50, main = "Tau Distribution", xlab = "Tau")
hist(spm_values, breaks = 50, main = "SPM Distribution", xlab = "SPM")
```

### 6.3 Tissue-Specificity Conservation

```r
# Classify genes by tissue-specificity
classify_specificity <- function(tau) {
    if (is.na(tau)) return("NA")
    if (tau < 0.3) return("Ubiquitous")
    if (tau < 0.7) return("Moderate")
    return("Specific")
}

tau_comparison$specificity_human <- sapply(tau_comparison$tau_human, classify_specificity)
tau_comparison$specificity_mouse <- sapply(tau_comparison$tau_mouse, classify_specificity)

# Cross-tabulation
table(tau_comparison$specificity_human, tau_comparison$specificity_mouse)

# Identify conservation categories
tau_comparison$specificity_conserved <- 
    tau_comparison$specificity_human == tau_comparison$specificity_mouse
```

---

## Phase 7: Expression Conservation Analysis

### 7.1 Pearson Correlation

```r
# Calculate expression correlation for ortholog pairs
# Using tissue mean expression

expression_correlation <- sapply(ortholog_gene_ids, function(gene) {
    human_profile <- tissue_mean_human[gene, ]
    mouse_profile <- tissue_mean_mouse[gene, ]
    
    if (any(is.na(human_profile)) || any(is.na(mouse_profile))) return(NA)
    
    cor(human_profile, mouse_profile, method = "pearson")
})

# Distribution of expression correlations
hist(expression_correlation, breaks = 50,
     main = "Expression Conservation (Pearson r)",
     xlab = "Correlation coefficient")

# Identify conserved vs divergent genes
conserved_genes <- names(which(expression_correlation > 0.8))
divergent_genes <- names(which(expression_correlation < 0.5))
```

### 7.2 Spearman Rank Correlation

```r
# Spearman correlation (non-parametric)
spearman_correlation <- sapply(ortholog_gene_ids, function(gene) {
    human_profile <- tissue_mean_human[gene, ]
    mouse_profile <- tissue_mean_mouse[gene, ]
    
    cor(human_profile, mouse_profile, method = "spearman")
})
```

### 7.3 Evolutionary Rate Estimation

```r
library(ape)
library(phytools)

# Load species tree with divergence times
species_tree <- read.tree("species_tree.nwk")

# Phylogenetic independent contrasts
# Requires divergence time-calibrated tree

# Calculate expression contrasts for each gene
pic_expression <- function(gene_expression, tree, species_names) {
    names(gene_expression) <- species_names
    pic(gene_expression, tree)
}

# Apply to all genes
expression_contrasts <- lapply(ortholog_gene_ids, function(gene) {
    gene_exp <- tissue_mean[gene, ]
    pic_expression(gene_exp, species_tree, species_names)
})
```

---

## Phase 8: Cross-Species Differential Expression

### 8.1 limma Analysis Framework

```r
library(limma)

# Create design matrix with species and tissue factors
design <- model.matrix(~ 0 + species + tissue + species:tissue)
colnames(design) <- make.names(colnames(design))

# Fit linear model
fit <- lmFit(corrected_exp, design)

# Define contrasts
# Example: tissue effect across species
contrast_matrix <- makeContrasts(
    tissueA_vs_tissueB = tissueA - tissueB,
    species_interaction = speciesspeciesA.tissueA - speciesspeciesB.tissueA,
    levels = design
)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract DE genes
de_genes <- topTable(fit2, coef = "tissueA_vs_tissueB", 
                     number = Inf, adjust.method = "BH")
```

### 8.2 DESeq2 for Count Data

```r
library(DESeq2)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_metadata,
    design = ~ species + tissue + species:tissue
)

# Run DESeq2
dds <- DESeq(dds)

# Extract results
# Tissue effect controlling for species
res <- results(dds, name = "tissue_tissueA_vs_tissueB")
```

---

## Phase 9: Visualization

### 9.1 Expression Heatmap

```r
library(pheatmap)

# Select variable genes
variable_genes <- names(sort(apply(corrected_exp, 1, var), decreasing = TRUE))[1:500]

# Z-score normalization
z_score <- t(scale(t(corrected_exp[variable_genes, ])))

# Heatmap
pheatmap(z_score,
         annotation_col = data.frame(species = batch, tissue = tissue),
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none")
```

### 9.2 Correlation Scatter Plot

```r
# Scatter plot of ortholog expression
ggplot(tau_comparison, aes(x = tau_human, y = tau_mouse)) +
    geom_point(alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(x = "Tau (Human)", y = "Tau (Mouse)",
         title = "Tissue-Specificity Conservation") +
    theme_bw()
```

### 9.3 Specificity Bar Plot

```r
# Top tissue-specific genes
top_specific <- head(sort(tau_values, decreasing = TRUE), 20)

# Bar plot
barplot(top_specific,
        main = "Top 20 Tissue-Specific Genes",
        ylab = "Tau Index",
        col = "steelblue",
        las = 2,
        cex.names = 0.7)
```

---

## FAQ

### Q1: No clear clustering by tissue after batch correction

**Possible causes:**
- Tissue samples are too different across species
- Batch correction overcorrected
- Tissue annotation mismatch

**Solutions:**
- Verify tissue homology between species
- Try different batch correction methods (ComBat vs MNN)
- Use less aggressive correction parameters

### Q2: Very low expression correlation between orthologs

**Possible causes:**
- Incorrect ortholog mapping
- Large evolutionary distance
- Different tissue sampling

**Solutions:**
- Re-verify ortholog mapping with multiple methods
- Focus on conserved genes (housekeeping) as positive control
- Use rank-based correlation (Spearman) instead

### Q3: ComBat fails with error

**Common fixes:**
- Ensure >3 samples per batch (species)
- Remove genes with zero variance
- Check for NA values in expression matrix

---

## Output Checklist

- [ ] Ortholog expression matrix (corrected)
- [ ] Tissue-specificity scores (tau/SPM)
- [ ] Expression correlation values
- [ ] PCA/t-SNE/UMAP visualizations
- [ ] Heatmap of variable genes
- [ ] Scatter plots of specificity conservation
- [ ] Summary statistics table
