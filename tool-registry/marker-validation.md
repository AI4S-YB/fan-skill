# Marker Validation

**Goal:** Verify marker specificity and polymorphism before wet-lab testing
**Best for:** All designed markers before ordering primers

## In-Silico PCR

```bash
# Check if primers amplify a unique product
blastn -query primers.fa -db reference_genome -outfmt 6 -task blastn-short
# Expect exactly 2 hits (F + R) within expected product size range
```

## Polymorphism Check

```r
# If you have population genotype data
check_polymorphism <- function(marker_genotypes) {
  maf <- mean(marker_genotypes) / 2
  missing_rate <- mean(is.na(marker_genotypes))
  list(maf = maf, missing_rate = missing_rate, usable = maf > 0.05 & missing_rate < 0.2)
}
```

## Plant-Specific Notes
- Polyploids: validate that primers amplify a SINGLE locus (not all homeologs)
- For breeding markers: MAF should be > 0.1 in the target breeding population
