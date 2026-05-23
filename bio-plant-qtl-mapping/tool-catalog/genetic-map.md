# Genetic Map Construction — est.map, ASMap, LepMap3

**Goal:** Build accurate genetic linkage maps from marker segregation data
**Best for:** All biparental populations; method selection depends on marker count and population type

## Prerequisites
- Clean genotype data (missing rate <20%, segregation distortion <20%)
- Population type information (F2, RIL, DH, BC)
- Chromosome/group assignments for each marker

## Method Selection Guide

| Method | Best For | Marker Limit | Speed |
|--------|----------|-------------|-------|
| R/qtl est.map | Standard populations, <5000 markers | ~5000 | Moderate |
| ASMap (R) | MAGIC, complex populations | ~50000 | Fast |
| LepMap3 | Ultra-dense maps (>10000 markers) | >100000 | Fast (Java) |

## Method 1: R/qtl est.map

```r
library(qtl)

# Read cross data
mycross <- read.cross("csv", file = "genotype.csv",
                      genotypes = c("AA", "AB", "BB"))

# Remove distorted and duplicate markers
mycross <- drop.nullmarkers(mycross)
gt <- geno.table(mycross)
bad_markers <- rownames(gt)[gt$P.value < 0.01]  # Segregation distortion p<0.01
mycross <- drop.markers(mycross, bad_markers)

# Build map
newmap <- est.map(mycross, error.prob = 0.001,
                  map.function = "kosambi",
                  maxit = 10000, tol = 1e-6)

# Replace old map with new map
mycross <- replace.map(mycross, newmap)

# Check map quality
plotMap(mycross, show.marker.names = FALSE)

# Check chromosome lengths
chrlen <- sapply(chrlen(mycross), function(x) max(x) - min(x))
chrlen
```

### Key R/qtl Parameters

| Parameter | Description | Recommended |
|-----------|-------------|-------------|
| error.prob | Genotyping error rate | 0.001 (SSR), 0.01 (GBS) |
| map.function | Map function | "kosambi" (plants), "haldane" (no interference) |
| maxit | Maximum EM iterations | 10000 (default), increase if warning |
| tol | Convergence tolerance | 1e-6 (default) |

## Method 2: ASMap (MAGIC / Complex Populations)

```r
library(ASMap)

# Convert R/qtl cross to ASMap mstmap object
msmap <- mstmap(mycross, bychr = TRUE, dist.fun = "kosambi",
                p.value = 1e-10,  # Linkage grouping threshold
                noMap.dist = 15,  # Max distance for clustering
                noMap.size = 0)

# Pull map
asmap <- pull.map(msmap)

# Compare with original map
compare_maps <- compare.order(mycross, msmap)

# Diagnostic: heatmap of pairwise recombination fractions
heatMap(msmap, chr = "all", lmax = 70)
```

### When to Use ASMap

- **MAGIC/NAM populations**: ASMap handles multi-parental populations robustly
- **Large marker sets**: Faster than est.map for >2000 markers
- **Segregation distortion**: Better tolerance of distorted markers

## Method 3: LepMap3 (Ultra-Dense Maps)

```bash
# Step 1: Prepare input — pedigree file
cat > pedigree.txt << 'EOF'
# family line1 line2
CHR POS P1 P2
fam1 1 2
EOF

# Step 2: Filter markers by segregation distortion
java -cp lepmap3.jar Filtering2 \
  data=genotype_called.txt \
  dataTolerance=0.001

# Step 3: Separate chromosomes (linkage groups)
java -cp lepmap3.jar SeparateChromosomes2 \
  data=genotype_called.txt \
  lodLimit=10 \
  sizeLimit=5

# Step 4: Order markers within each linkage group
java -cp lepmap3.jar OrderMarkers2 \
  map=map.txt \
  data=genotype_called.txt \
  useMorgan=1 \
  improveOrder=1

# Step 5: Join maps from individual linkage groups
java -cp lepmap3.jar JoinSingles2All \
  map=map.txt \
  data=genotype_called.txt \
  lodLimit=5 \
  iterate=1
```

### When to Use LepMap3

- **>10000 markers**: R/qtl est.map becomes too slow
- **GBS/WGS data**: Handles massive missing data patterns well
- **Complex pedigrees**: Supports multi-family designs
- **Sex-averaged and sex-specific maps**: Can model sex-specific recombination

## Map Quality Checks

```r
# Check marker order consistency
plotRF(mycross)  # Should show decreasing RF with distance

# Check chromosome lengths
# For most plant species, genetic maps are 1000-3000 cM total
total_map_length <- sum(chrlen(mycross))

# Check marker spacing
marker_density <- nmarkers(mycross) / chrlen(mycross)

# Identify potential mis-ordered markers
# Large gaps (>30 cM) may indicate missing data or ordering errors
```

## Plant-Specific Notes

- **Map expansion**: Plant genetic maps from R/qtl est.map are often 10-30% longer than expected due to genotyping errors. This is normal.
- **Kosambi vs Haldane**: Use Kosambi for plants — it accounts for crossover interference, which is present in most plant species.
- **Centromere location**: Look for regions with suppressed recombination (tightly clustered markers at the genetic level but spaced widely at the physical level).
- **Polyploidy**: For polyploid species, each subgenome may have different map lengths. Allopolyploids (e.g., wheat) need separate maps per subgenome.
- **GBS data**: Expect 30-70% missing data. LepMap3 handles this better than R/qtl.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `"est.map: markers at the same position"` | Duplicate markers after error correction | Use `jittermap()` to add small jitter |
| `"est.map: negative distances"` | Markers with conflicting recombination information | Check `plotRF()` for unusually high RF values |
| `"est.map did not converge"` | Too many markers or genotyping errors | Increase maxit or reduce error.prob |
| Map length >5000 cM | Massive genotyping errors or wrong marker order | Check segregation distortion; remove problematic markers |
| Chromosome fragments (multiple small LGs for one chr) | LOD threshold too high in grouping | Lower threshold or manually merge groups |
