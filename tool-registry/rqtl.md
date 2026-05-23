# R/qtl Core — Data Import, Map Visualization, Single-QTL Scan

**Goal:** Import genotype+phenotype data, visualize genetic maps, run single-QTL interval mapping
**Best for:** All standard biparental populations (F2, RIL, DH, BC)

## Prerequisites
- R 4.0+, qtl package (1.50+)
- Genotype data: CSV format (markers x individuals)
- Phenotype data: CSV format (individuals x traits)
- Population information: type, number of chromosomes, marker positions

## Basic Usage

### Data Import

```r
library(qtl)

# Read cross data — format depends on input files
# Option A: CSV format (most common for plant data)
mycross <- read.cross(
  format = "csv",
  file = "genotype.csv",
  phefile = "phenotype.csv",
  genotypes = c("AA", "AB", "BB"),  # F2 coding
  na.strings = c("-", "NA")
)

# Option B: Mapmaker format (legacy plant datasets)
# mycross <- read.cross(format = "mm", file = "data.raw")

# Option C: R/qtl format (RDS saved cross object)
# mycross <- readRDS("cross.rds")

# Summary of the cross
summary(mycross)

# Plot missing genotype pattern
plotMissing(mycross)
```

### Population Type Handling

```r
# Convert to appropriate cross type
class(mycross)  # Should be one of: "f2", "bc", "riself", "dh"

# F2 → RIL conversion (if selfing advanced)
# mycross <- convert2riself(mycross)

# Verify marker coding
plotGeno(mycross, chr = 1, ind = 1:20)
```

### Genetic Map Visualization

```r
# Plot full genetic map
plotMap(mycross, main = "Genetic Map")

# Plot single chromosome
plotMap(mycross, chr = "1", show.marker.names = TRUE)

# Check map quality
plotRF(mycross)  # Recombination fraction heatmap
```

### Interval Mapping (IM)

```r
# Calculate genotype probabilities
mycross <- calc.genoprob(mycross, step = 1, error.prob = 0.001)

# Single-QTL scan via interval mapping
out.im <- scanone(mycross, method = "em", pheno.col = "trait1")

# View results
summary(out.im)
plot(out.im)

# Permutation test for LOD threshold
operm <- scanone(mycross, method = "em", pheno.col = "trait1",
                 n.perm = 1000)
summary(operm)
lod_threshold <- summary(operm)[1]  # 5% genome-wide threshold

# Identify significant QTL
summary(out.im, perms = operm, alpha = 0.05, pvalues = TRUE)
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| step | 1 cM (dense map) or 2 cM (sparse map) | Step size for genotype probability calculation |
| error.prob | 0.001 for SSR, 0.01 for GBS | Genotyping error rate assumption |
| method | "em" (standard) or "hk" (Haley-Knott, faster) | EM more accurate, HK faster for large datasets |
| n.perm | 1000 for n>=100, 500 for n<100 | Permutation count for LOD threshold |

## Plant-Specific Notes

- **F2 coding**: Plant F2 populations use AA/AB/BB coding. Make sure `genotypes` parameter matches your data — some plant data use 0/1/2 or A/H/B.
- **Dominance effects**: F2 populations can estimate both additive and dominance effects via `scanone(..., model = "binary")` or `sim.geno`.
- **RIL selfing**: For RILs derived by selfing (not sib-mating), use `class(mycross) <- "riself"`. If derived by sib-mating, use `"risib"`.
- **Segregation distortion**: Run `geno.table(mycross)` to check for markers deviating from expected Mendelian ratios.
- **Map expansion**: Genetic maps from plant populations are often longer than expected due to genotyping errors or high recombination in certain regions.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `"Inconsistent number of individuals"` | Mismatch between genotype and phenotype | Check sample IDs in both files |
| `"Missing chromosomes"` | Marker data missing chromosome information | Add chromosome column to input CSV |
| `"Genotype probabilities not calculated"` | Forgot calc.genoprob | Run `calc.genoprob()` before scanone |
| `"Error in est.map: negative distances"` | Markers with identical or conflicting positions | Use `drop.nullmarkers()` and `drop.dupmarkers()` |
| `"No significant QTL found"` | Low power or trait has no major QTL | Report as null result; lower alpha to suggestive threshold |
