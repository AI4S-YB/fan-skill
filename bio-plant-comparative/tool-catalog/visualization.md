# Comparative Genomics Visualization

## Synteny Dot Plot

```r
library(ggplot2)
ggplot(synteny_blocks, aes(x = speciesA_pos, y = speciesB_pos, color = chromosome)) +
  geom_point(size = 0.1, alpha = 0.5) +
  facet_grid(chrB ~ chrA, scales = "free") +
  labs(x = "Species A", y = "Species B") + theme_minimal()
```

## Ks Distribution with Peak Fitting

```r
ggplot(ks_data, aes(x = Ks)) +
  geom_histogram(bins = 100, fill = "steelblue") +
  geom_vline(xintercept = peak_positions, color = "red", linetype = "dashed") +
  labs(x = "Ks (synonymous substitutions per site)", y = "Gene pair count")
```

## dN/dS Distribution

```r
ggplot(dnds_data, aes(x = dNdS)) +
  geom_density(fill = "grey", alpha = 0.5) +
  geom_vline(xintercept = 1, color = "red") +
  labs(x = "dN/dS ratio", y = "Density")
```
