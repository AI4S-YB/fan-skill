# DADA2 Denoising

**Goal:** Infer amplicon sequence variants (ASVs) from raw sequencing reads with single-nucleotide resolution
**Best for:** All Illumina amplicon data — recommended over OTU clustering for modern sequencing

## Prerequisites
- R 4.0+
- R packages: dada2, ShortRead, Biostrings, ggplot2
- QIIME 2 (optional, for integration with downstream analyses)

## Standard DADA2 Workflow in R

```r
library(dada2)

# Define paths
path <- "raw_reads/"
fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# Quality profiles
plotQualityProfile(fnFs[1:2])  # Check forward reads
plotQualityProfile(fnRs[1:2])  # Check reverse reads

# Filter and trim
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(240, 200),
                     maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = 16)

# Learn error rates
errF <- learnErrors(filtFs, multithread = 16)
errR <- learnErrors(filtRs, multithread = 16)
plotErrors(errF, nominalQ = TRUE)

# Sample inference (dereplication + denoising)
dadaFs <- dada(filtFs, err = errF, multithread = 16, pool = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = 16, pool = TRUE)

# Merge paired-end reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                      minOverlap = 12, maxMismatch = 0)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                     multithread = 16, verbose = TRUE)

# Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names
print(track)
```

## Key Parameters

| Parameter | Purpose | Typical Value |
|-----------|---------|---------------|
| truncLen | Trim reads to fixed length | Based on quality plot |
| maxEE | Maximum expected errors per read | 2 for both directions |
| truncQ | Truncate at first base with Q < value | 2 |
| maxN | Maximum ambiguous bases | 0 |
| minOverlap | Minimum overlap for merging | 12 |
| pool | Pool samples for sample inference | TRUE (better for rare variants) |
| rm.phix | Remove PhiX spike-in reads | TRUE |

## QIIME 2 DADA2 Plugin

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads 16 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza
```

## Plant-Specific Considerations

- Plant chloroplast and mitochondrial 16S reads will be abundant in root/leaf samples — do NOT filter them during DADA2; remove them after taxonomy assignment
- If using PNA clamps, expect lower total reads but higher bacterial proportion
- Rhizosphere samples typically have high diversity; `pool = TRUE` helps detect rare ASVs
- For low-biomass samples (endosphere, phyllosphere), be conservative with `maxEE` to avoid filtering rare real sequences

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No reads passed the filter" | truncLen too long or quality too low | Check quality profiles; reduce truncLen |
| "Error in dada: $err" | Too few reads for error learning | Use pool=TRUE or increase sample size |
| "All samples have zero sequences" | Mismatched read pairs or wrong file paths | Verify file naming pattern |
| Memory exhausted | Too many samples loaded at once | Process in batches or increase RAM |
| "mergePairs: no pairs merged" | Overlap too short or read orientation wrong | Check if reads are in correct orientation |
