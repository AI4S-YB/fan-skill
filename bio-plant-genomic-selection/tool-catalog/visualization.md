# Genomic Selection Visualization

**Goal:** Accuracy boxplot, predicted vs observed, GEBV distribution, selection differential
**Best for:** All GS results

## Accuracy Boxplot

```r
library(ggplot2)

ggplot(cv_results, aes(x = method, y = accuracy)) +
  geom_boxplot(fill = "steelblue") +
  geom_jitter(width = 0.1, alpha = 0.3) +
  labs(y = "Prediction Accuracy (Pearson r)", x = "") +
  theme_minimal()
```

## Predicted vs Observed

```r
ggplot(predictions, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = min(observed), y = max(predicted),
           label = paste0("r = ", round(accuracy, 3))) +
  labs(x = "Observed Phenotype", y = "GEBV") +
  theme_minimal()
```

## GEBV Distribution with Selection

```r
# Highlight top 10% selected individuals
top_cutoff <- quantile(gebvs, 0.90)
ggplot(data.frame(gebv = gebvs), aes(x = gebv)) +
  geom_histogram(aes(fill = gebv >= top_cutoff), bins = 30) +
  scale_fill_manual(values = c("grey", "red"), labels = c("Not selected", "Selected")) +
  labs(x = "GEBV", fill = "") +
  theme_minimal()
```

## Plant-Specific Notes

- Multi-trait GS: add scatter plot of trait2 vs trait1 GEBVs, colored by selection status
- Report selection differential (mean GEBV of selected − mean GEBV of all)
