# Cross-Validation for Genomic Selection

**Goal:** Estimate prediction accuracy without an independent validation set
**Best for:** All GS analyses — required before applying GS in breeding

## Standard k-Fold CV

```r
cv_gs <- function(geno, pheno, k = 5, repeats = 5) {
  n <- length(pheno)
  accuracies <- c()

  for (rep in 1:repeats) {
    folds <- sample(rep(1:k, length.out = n))
    for (i in 1:k) {
      test_idx <- which(folds == i)
      train_idx <- which(folds != i)

      # Train on training set
      model <- mixed.solve(y = pheno[train_idx], K = G[train_idx, train_idx])

      # Predict test set
      pred <- model$u[test_idx]
      accuracies <- c(accuracies, cor(pred, pheno[test_idx]))
    }
  }
  return(mean(accuracies))
}
```

## Stratified CV for Structured Populations

```r
# Don't use random folds when you have subpopulations!
# Instead: leave-one-subpopulation-out
for (pop in unique(population_labels)) {
  test_idx <- which(population_labels == pop)
  train_idx <- which(population_labels != pop)
  # ... train and predict
}
```

## Plant-Specific CV Design

| Scenario | CV Strategy |
|----------|------------|
| Unrelated individuals | 5-fold × 5 repeats |
| Full-sib families | Leave-one-family-out |
| Multi-environment | Leave-one-environment-out |
| Breeding pipeline (generations) | Forward prediction (old → new) |
