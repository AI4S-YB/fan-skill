# Marker Visualization

## Primer Position Schematic

```r
library(ggplot2)
ggplot(primer_data, aes(x = pos, y = 1)) +
  geom_segment(aes(xend = pos + length, yend = 1, color = type), size = 3) +
  geom_vline(xintercept = snp_pos, linetype = "dashed", color = "red") +
  labs(x = "Genomic position (bp)", y = "") + theme_minimal()
```

## Marker Polymorphism Plot

```r
ggplot(pop_data, aes(x = allele, fill = allele)) +
  geom_bar() + facet_wrap(~ marker) +
  labs(title = "Marker Polymorphism in Breeding Population")
```

## Parental Combination Matrix

```r
pheatmap::pheatmap(compatibility_matrix,
  main = "Parental Pair Complementary Scores",
  color = colorRampPalette(c("white", "darkgreen"))(100))
```
