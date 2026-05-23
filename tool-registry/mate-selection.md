# Mate Selection (Optimal Contribution Selection)

**Goal:** Select optimal parental combinations for the next breeding cycle with genetic gain and diversity constraints
**Best for:** Designing crossing blocks in plant hybrid breeding programs

## Prerequisites

- R 4.0+ with optiSel, AlphaMate
- Estimated breeding values (EBVs) for all candidate parents
- Pedigree or genomic relationship matrix (A or G matrix)
- Breeding constraints (budget, number of crosses, diversity targets)

## Optimal Contribution Selection

```r
library(optiSel)

# Prepare candidate data
candidates <- data.frame(
  Indiv = inbred_lines$ID,
  Sex = inbred_lines$Sex,  # "male" or "female" for hybrid breeding
  EBV = inbred_lines$Yield_EBV,
  Group = inbred_lines$heterotic_group
)

# Build pedigree object
ped <- prePed(candidates, Pedigree)

# Calculate coancestry
# For hybrid breeding: males from group A, females from group B
males <- candidates[candidates$Sex == "male", ]
females <- candidates[candidates$Sex == "female", ]

# OCS for male parents
cand_male <- candes(phen = males, ped = ped)
con_male <- list(
  uniform = "female",  # Equal contributions to maintain size
  ub = rep(0.05, nrow(males))  # Max 5% contribution per male
)

# Solve OCS problem
opt_male <- opticont(
  method = "min.IBD",
  cand = cand_male,
  con = con_male,
  trace = FALSE
)

summary(opt_male)
```

## Genetic Algorithm for Mate Allocation

```r
# Define fitness function
# Maximize: genetic gain - lambda * inbreeding
mate_fitness <- function(mating_plan, ebv, kinship, lambda = 0.5) {
  n_crosses <- nrow(mating_plan)
  genetic_gain <- mean(ebv[mating_plan$Male] + ebv[mating_plan$Female]) / 2
  inbreeding <- mean(kinship[mating_plan$Male, mating_plan$Female])
  return(genetic_gain - lambda * inbreeding)
}

# Simple random search with constraints
set.seed(42)
best_plan <- NULL
best_score <- -Inf

for (iter in 1:10000) {
  # Randomly select male-female pairs
  male_sample <- sample(1:length(male_ids), n_crosses, replace = TRUE)
  female_sample <- sample(1:length(female_ids), n_crosses, replace = TRUE)

  # Constraint: no self-crosses (different groups)
  mating_plan <- data.frame(
    Male = male_ids[male_sample],
    Female = female_ids[female_sample]
  )

  # Constraint: max 3 crosses per parent
  male_counts <- table(mating_plan$Male)
  female_counts <- table(mating_plan$Female)
  if (any(male_counts > 3) || any(female_counts > 3)) next

  # Compute fitness
  score <- mate_fitness(mating_plan, ebvs_hybrid, kinship)
  if (score > best_score) {
    best_score <- score
    best_plan <- mating_plan
  }
}

print(paste("Best fitness score:", round(best_score, 3)))
print(paste("Selected", nrow(best_plan), "crosses"))
```

## Breeding Constraints

```r
# Practical constraints for plant breeding
constraints <- list(
  max_crosses_per_male = 5,     # Pollen limitation
  max_crosses_per_female = 3,   # Seed production capacity
  min_total_parents = 20,       # Maintain genetic diversity
  max_relatedness = 0.125,      # Avoid half-sib matings (F < 0.125)
  heterotic_group_rule = "between",  # Only cross between groups
  elite_reserve = 2,            # Reserve 2 top parents for backup
  total_crosses = 100           # Budget constraint
)

# Apply constraints step-by-step
plan_filtered <- all_possible_crosses

# 1. Heterotic group constraint
plan_filtered <- plan_filtered[
  plan_filtered$Male_Group != plan_filtered$Female_Group,
]

# 2. Relatedness constraint
plan_filtered <- plan_filtered[
  kinship_cross < 0.125,
]

# 3. Select based on predicted performance
plan_filtered$predicted_yield <- predict_hybrid(plan_filtered)
plan_filtered <- plan_filtered[order(-plan_filtered$predicted_yield), ]

# 4. Apply per-parent caps
selected_crosses <- data.frame()
male_used <- c()
female_used <- c()

for (i in 1:nrow(plan_filtered)) {
  cross <- plan_filtered[i, ]
  if (male_used[cross$Male] %||% 0 >= constraints$max_crosses_per_male) next
  if (female_used[cross$Female] %||% 0 >= constraints$max_crosses_per_female) next

  selected_crosses <- rbind(selected_crosses, cross)
  male_used[cross$Male] <- (male_used[cross$Male] %||% 0) + 1
  female_used[cross$Female] <- (female_used[cross$Female] %||% 0) + 1

  if (nrow(selected_crosses) >= constraints$total_crosses) break
}
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| Max crosses per parent | 3-5 | Practical pollen/seed constraints |
| Lambda (gain vs inbreeding) | 0.3-0.7 | Lower = prioritize gain, higher = prioritize diversity |
| Min parent representation | >= 10 males, >= 10 females | Avoid narrowing genetic base |
| Reserve parents | 2-3 per group | Backup for failed crosses |
| Target inbreeding rate | DeltaF < 0.01 per generation | Sustainable long-term breeding |

## Plant-Specific Notes

- For hybrid crops (maize, rice, sorghum): select distinct male and female pools from heterotic groups
- For synthetic varieties (alfalfa): balance parent number with gain
- Recurrent selection programs: use OCS iteratively each cycle
- For genomic selection-integrated programs: update EBVs and re-optimize annually
- The "genetic gain - diversity" tradeoff is fundamental: accepting ~5% less gain per cycle can double the effective population size

## Genetic Gain Projection

```r
# Project genetic gain over cycles
n_cycles <- 10
gain_per_cycle <- mean(selected_ebvs) - mean(population_ebvs)
cumulative_gain <- cumsum(rep(gain_per_cycle, n_cycles))

plot(1:n_cycles, cumulative_gain, type = "b",
     xlab = "Breeding Cycle", ylab = "Cumulative Genetic Gain",
     main = "Projected Genetic Gain")
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Selected parents all related | Lambda too low (prioritized gain) | Increase lambda, check kinship matrix |
| No feasible solution | Too many constraints | Relax per-parent caps or min parent requirement |
| Genetic gain negative | Model predicts poorly with these parents | Validate predictions on test data |
| All top crosses use same parents | Distribution of EBVs very skewed | Add diversity constraints or weighted selection |
