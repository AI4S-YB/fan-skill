# Spatial Correction for Field Trials

**Goal:** Correct for spatial trends in field trials using row-column models
**Best for:** Field trials with known row/column positions, CV > 15%

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| Model type | SpATS (splines) | Regular grid layout with few rows/columns: AR1xAR1; irregular spacing or very large trials (>5000 plots): SpATS | AR1xAR1 requires separable row and column autocorrelation; SpATS handles irregular spatial trends and field gradient complexity more flexibly |
| nseg (SpATS) | c(10, 10) | Small trials (<100 plots): c(5, 5); large trials (>1000 plots): c(20, 20) | Segment number balances flexibility vs overfitting; too many segments capture genetic signal as spatial noise and reduce genotype variance |
| Residual check | Visual variogram | Formal hypothesis test needed for publication: use Moran's I on model residuals | Significant residual spatial autocorrelation indicates underspecified spatial model; switch from AR1 to SpATS or add nugget effect |
| Trial design recognition | Lattice / alpha design | RCBD with no spatial blocking: always add spatial correction; augmented design with sparse checks: use SpATS with fixed check term | Incomplete block designs already partially capture spatial variation; spatial correction recovers additional precision for genotype BLUPs |

## Prerequisites
- R 4.0+, sommer or SpATS
- Row and column coordinates for each plot

## AR1×AR1 Spatial Model (sommer)

```r
library(sommer)
ans <- mmer(trait ~ genotype,
  random = ~ vs(row) + vs(col) + vs(units),
  rcov = ~ vs(units),
  data = pheno)
```

## SpATS (Spatial Analysis of Field Trials with Splines)

```r
library(SpATS)
ans <- SpATS(response = "trait", genotype = "genotype",
  spatial = ~ PSANOVA(row, col, nseg = c(10, 10)),
  data = pheno)
```

## Plant-Specific Notes
- Spatial variation is common in large breeding trials — always check CV
- Row-column designs (lattice, alpha) already partially control for spatial effects
- Spatial correction can substantially increase heritability estimates

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "Model did not converge" | Too few rows/columns or extreme imbalance | Use spline-based SpATS instead |
| Heritability estimate decreased after correction | Spatial correction captured genetic signal as noise | Verify genotype randomization in field; check that spatial terms do not confound with genetic terms; reduce nseg |
| Residuals still show spatial pattern after correction | Model underspecified or non-stationary spatial trend | Try adding nugget effect to AR1 model; use anisotropic SpATS with column-specific nseg; check for edge effects |
