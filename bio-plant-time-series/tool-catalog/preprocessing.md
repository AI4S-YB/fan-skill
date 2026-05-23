# Preprocessing: Expression Smoothing

**Goal:** Denoise time-series expression data while preserving true temporal
signal (peaks, inflection points, transition timing).

**Best for:** Normalized expression matrices (DESeq2 VST, rlog, or
variance-stabilized) with 3+ time points.

## Prerequisites

- R 4.0+, `stats` (base R), `psych` or `Mfuzz` filtering functions
- Input: gene x time-point matrix, genes in rows, time points in columns
- Normalized counts (not raw) -- DESeq2 VST or rlog recommended
- Metadata: vector of time-point labels (numeric, in order)

## Spline Smoothing (>= 5 time points)

```r
# Fit smoothing spline per gene with GCV df selection
smooth_gene_spline <- function(expr_vec, time_points) {
  fit <- smooth.spline(time_points, expr_vec, cv = TRUE)
  return(predict(fit, time_points)$y)
}

# Apply to whole matrix
smoothed_mat <- t(apply(expr_matrix, 1, function(row) {
  smooth_gene_spline(row, time_points)
}))

# Check autocorrelation of residuals
residuals <- expr_matrix - smoothed_mat
acf_vals <- apply(residuals, 1, function(r) {
  acf(r, plot = FALSE, lag.max = 1)$acf[2]
})
cat("Median residual ACF(1):", median(abs(acf_vals), na.rm = TRUE), "\n")
# Acceptable: < 0.3
```

## LOESS Smoothing (< 5 time points)

```r
# LOESS per gene with small span for sparse data
smooth_gene_loess <- function(expr_vec, time_points, span = 0.75) {
  if (length(time_points) < 4) span <- 1.0  # max smooth for very few points
  fit <- loess(expr_vec ~ time_points, span = span)
  return(predict(fit, time_points))
}

smoothed_mat <- t(apply(expr_matrix, 1, function(row) {
  smooth_gene_loess(row, time_points)
}))
```

## Parameters

| Parameter | Spline default | LOESS default | Description |
|-----------|---------------|---------------|-------------|
| df / spar | GCV auto | -- | Degrees of freedom (spline); NULL = cross-validation |
| span | -- | 0.75 | LOESS smoothing span; 0-1, higher = smoother |
| cv | TRUE | -- | Use generalized cross-validation for df selection |
| boundary | natural | -- | Boundary condition for spline |

## Plant-Specific Notes

- **Circadian data**: do NOT smooth before periodicity detection.
  Smooth after confirming a gene is non-periodic, or use a periodic spline.
- **Developmental series**: spline df of n_timepoints/2 often works well.
  Do not exceed n_timepoints - 1.
- **Stress time courses**: rapid induction (0.5-2 hr) can be smoothed
  away by aggressive LOESS. Use span >= 0.5 to preserve sharp peaks.
- **Diurnal field data**: temperature fluctuation can introduce
  24-hr oscillation in smoothing residuals. Check before proceeding.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| `smooth.spline: NA/NaN values` | Missing expression values | Impute missing time points or remove gene |
| `loess: span too small` | Fewer data points than span allows | Increase span or use `simple` surface |
| Residual autocorrelation > 0.5 | Over-smoothing removed real signal | Decrease df (spline) or span (LOESS) |
| Boundary values diverge | Spline boundary condition too loose | Set `all.knots = TRUE` or use natural spline |
| Zero variance after smoothing | Gene has constant expression | Remove constant genes before clustering |
