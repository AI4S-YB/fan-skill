# Phenotype Analysis Visualization

## Heritability Bar Plot
```r
ggplot(h2_data, aes(x = trait, y = H2)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = H2 - se, ymax = H2 + se), width = 0.2) +
  labs(y = "Broad-sense Heritability (H²)", x = "") + theme_minimal()
```

## BLUP Ranking Plot
```r
blup_df <- blup_df[order(blup_df$blup), ]
blup_df$rank <- 1:nrow(blup_df)
ggplot(blup_df, aes(x = rank, y = blup)) +
  geom_point() + geom_errorbar(aes(ymin = blup - se, ymax = blup + se), alpha = 0.3) +
  labs(x = "Genotype Rank", y = "BLUP") + theme_minimal()
```

## GGE Biplot
Use `GGEBiplots` package for publication-quality biplots.

## Spatial Heatmap
```r
ggplot(pheno, aes(x = col, y = row, fill = trait)) +
  geom_tile() + scale_fill_viridis_c() +
  labs(title = "Spatial Distribution of Trait Values") + theme_minimal()
```
