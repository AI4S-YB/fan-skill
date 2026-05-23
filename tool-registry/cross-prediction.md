# Cross Prediction (Genomic BLUP / BayesR)

**Goal:** Predict the performance of unobserved hybrid crosses using genomic data
**Best for:** Screening thousands of potential hybrid combinations without field testing

## Prerequisites

- R 4.0+ with BGLR, rrBLUP, or sommer
- Genotype matrix (SNP: inbred lines x markers)
- Phenotype data from observed hybrid crosses
- Mating design information (which parent combinations were tested)

## GBLUP for Hybrid Prediction

```r
library(sommer)

# Step 1: Build genomic relationship matrices
# Additive G matrix
G_A <- A.mat(SNP_matrix)

# Step 2: Set up hybrid prediction model
# Model: y = mu + GCA_male + GCA_female + SCA

# Build design matrices for GCA effects
Z_male <- model.matrix(~ 0 + Male, data = hybrid_data)
Z_female <- model.matrix(~ 0 + Female, data = hybrid_data)

# Fit model
model_hybrid <- mmer(
  fixed = Yield ~ Env,
  random = ~ vsr(Male, Gu = G_A) +
             vsr(Female, Gu = G_A) +
             Male:Female,
  rcov = ~ units,
  data = hybrid_data
)

# Step 3: Predict untested crosses
# Create prediction dataset with all possible crosses
all_crosses <- expand.grid(
  Male = unique(inbred_lines$ID),
  Female = unique(inbred_lines$ID)
)
# Remove self-crosses and tested crosses
all_crosses <- all_crosses[all_crosses$Male != all_crosses$Female, ]
untested <- all_crosses[!paste(all_crosses$Male, all_crosses$Female) %in%
                         paste(hybrid_data$Male, hybrid_data$Female), ]

# Predict
predictions <- predict(model_hybrid, newdata = untested)
```

## BayesR for Hybrid Prediction

```r
library(BGLR)

# Prepare response
y <- hybrid_data$Yield

# Prepare incidence matrices for GCA effects
# For each inbred line, calculate hybrid genetic value
# G_hybrid = 0.5 * (g_parent1 + g_parent2) + d_cross

ETA <- list(
  list(X = Z_male, model = "BRR"),     # GCA male (Bayes Ridge Regression)
  list(X = Z_female, model = "BRR"),   # GCA female
  list(X = Z_cross, model = "BRR")     # SCA
)

# Fit model
fit_bayes <- BGLR(
  y = y,
  ETA = ETA,
  nIter = 20000,
  burnIn = 5000,
  thin = 5,
  saveAt = "bayesR_"
)

# Extract predictions
pred_train <- fit_bayes$yHat
```

## Cross Validation for Prediction Accuracy

```r
# K-fold cross-validation by cross (not by parent)
# This is critical: randomly split HYBRID COMBINATIONS, not parents

set.seed(42)
n_crosses <- nrow(hybrid_data)
folds <- sample(rep(1:5, length.out = n_crosses))

predicted <- numeric(n_crosses)

for (k in 1:5) {
  train_idx <- which(folds != k)
  test_idx <- which(folds == k)

  # Train on 4 folds
  model_cv <- mmer(
    fixed = Yield ~ 1,
    random = ~ vsr(Male, Gu = G_A) + vsr(Female, Gu = G_A),
    rcov = ~ units,
    data = hybrid_data[train_idx, ]
  )

  # Predict test fold
  predicted[test_idx] <- predict(model_cv,
                                  newdata = hybrid_data[test_idx, ])$predicted
}

# Accuracy
prediction_accuracy <- cor(hybrid_data$Yield, predicted, method = "pearson")
print(paste("Prediction accuracy (r):", round(prediction_accuracy, 3)))
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| Cross-validation folds | 5-10 | Balance bias-variance |
| Minimum crosses per parent | >= 5 | Stable GCA estimates for prediction |
| Minimum training crosses | >= 100 | Sufficient for genomic prediction |
| BayesR nIter | 20,000-50,000 | Adequate MCMC convergence |
| BayesR burnIn | 5,000-10,000 | Discard initial samples |

## Plant-Specific Notes

- Prediction accuracy is lower for traits with strong GxE interaction (e.g., yield)
- Prediction accuracy improves when using both additive and dominance effects
- For crops with few historical crosses: use pedigree-based prediction instead
- Rice: GCA-based prediction is often sufficient (r = 0.3-0.5)
- Maize: need dominance effects for accurate prediction (r = 0.5-0.7)
- Multi-environment data substantially improves prediction accuracy

## Forward Prediction (Temporal Validation)

```r
# More realistic than random CV:
# Train on old crosses, predict new crosses

# Split by year
train_years <- c(2018, 2019, 2020)
test_years <- c(2021, 2022)

train_idx <- hybrid_data$Year %in% train_years
test_idx <- hybrid_data$Year %in% test_years

# Train
model_fwd <- mmer(fixed = Yield ~ 1,
                  random = ~ vsr(Male, Gu = G_A) + vsr(Female, Gu = G_A),
                  rcov = ~ units,
                  data = hybrid_data[train_idx, ])

# Predict forward
pred <- predict(model_fwd, newdata = hybrid_data[test_idx, ])
forward_accuracy <- cor(hybrid_data$Yield[test_idx], pred$predicted)
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Prediction accuracy < 0 | Model worse than mean prediction | Check data quality, increase training size |
| Accuracy much higher in CV than forward | Data leakage (CV splits by parent instead of cross) | Split by HYBRID combination, not parent |
| Prediction of untested crosses all identical | Model missing SCA/dominance terms | Add specific combining ability random term |
| MCMC not converging | Insufficient iterations or poor mixing | Increase nIter, check trace plots |
