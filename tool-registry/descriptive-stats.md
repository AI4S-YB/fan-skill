# Descriptive Statistics for Phenotype Data

**Goal:** Summarize phenotype data with descriptive statistics and visualizations
**Best for:** Initial data exploration before any analysis

## Prerequisites
- R 4.0+, ggplot2, dplyr

## Basic Summary

```r
library(dplyr)
pheno %>%
  group_by(genotype) %>%
  summarise(
    mean = mean(trait, na.rm = TRUE),
    sd = sd(trait, na.rm = TRUE),
    cv = sd/mean * 100,
    n = n()
  )
```

## Distribution Plot

```r
library(ggplot2)
ggplot(pheno, aes(x = trait)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ environment) +
  labs(title = "Trait Distribution by Environment")
```

## Plant-Specific Notes
- Check for outliers before GWAS/GS — extreme values may be data entry errors
- Box-Cox transformation may be needed for non-normal traits
- Report CV as a measure of trial quality
