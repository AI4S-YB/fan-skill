# Parental Combination Recommendation

**Goal:** Recommend optimal parental combinations for crossing
**Best for:** Breeding decision support after GWAS/GS

## Complementary Allele Strategy (Pyramiding)

```r
# For each parental pair, count complementary favorable alleles
complementary_score <- function(parent1_alleles, parent2_alleles, favorable) {
  p1_good <- parent1_alleles == favorable
  p2_good <- parent2_alleles == favorable
  sum(p1_good | p2_good)  # Union of good alleles
}
```

## GEBV-Based Pairing

```r
# Cross parents with highest mid-parent GEBV
mid_parent_gebv <- (gebv[parent1] + gebv[parent2]) / 2
best_pairs <- order(mid_parent_gebv, decreasing = TRUE)[1:10]
```

## Plant-Specific Notes
- For self-pollinated crops: F1 -> F2 onward for selection
- For hybrid crops: test cross GCA before recommending
- Consider flowering time synchrony between parental lines
