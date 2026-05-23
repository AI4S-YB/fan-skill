# Plant Time-Series Expression: Analyst Notebook

This notebook helps you reason like an experienced plant
developmental biologist when analyzing time-series expression data.

## Before You Start: Understand Your Experimental Design

### What is the biological time scale?

Plant time-series experiments fall into distinct categories.
Knowing which one you are in shapes every downstream decision:

| Category | Time scale | Example |
|----------|-----------|---------|
| Diurnal / circadian | 2-4 hr intervals, 24-72 hr | Leaf sampling every 4 hr for 48 hr |
| Developmental gradient | Days to weeks | Leaf 1 through leaf 8 along the shoot |
| Stress response | Minutes to hours | 0, 0.5, 1, 2, 4, 8, 24 hr post-treatment |
| Tissue atlas | Stages / organs | Root, stem, leaf, flower, seed at multiple stages |
| Hormone response | 0.5-24 hr | Auxin or ABA time course |

Why this matters:
- **Circadian**: never smooth away the ~24 hr oscillation. Do not treat it
  as noise. Check periodicity with a periodogram BEFORE smoothing. A sine-wave
  pattern is biology, not artefact.
- **Developmental gradient**: expect monotonic or sigmoidal trends. Clusters
  dominated by "up then down" patterns may indicate sequential tissue sampling
  rather than true temporal regulation.
- **Stress response**: rapid induction peaks often happen at 2-4 hr. Make
  sure your time-point density captures the peak; otherwise you will miss
  the most informative window.

### How many time points do you have?

- **n < 5**: very coarse temporal resolution. You can detect large-effect
  patterns (constitutively on/off, monotonic up/down). Cluster counts
  should be small (k = 3-5). LOESS smoothing only. Be honest: you are
  seeing trends, not trajectories.
- **n = 5-10**: adequate for most developmental studies. Spline smoothing
  works. k = 6-12 for clustering. Can distinguish transient from sustained
  responses.
- **n > 10**: dense time course. Can use STEM for statistically rigorous
  trend detection. Spline df can be higher. You can resolve multi-phasic
  patterns (e.g., early up then late down).

### Tissue atlas interpretation

Plant tissue atlases (e.g., rice eFP browser, Arabidopsis TraVA, maize
Expression Atlas) provide spatial context. When interpreting clusters:

1. Cross-reference cluster members with known tissue markers.
2. A cluster enriched in "photosynthesis" GO terms that peaks at midday
   in leaf tissue is expected -- do not overinterpret.
3. Clusters where root-specific genes co-express with stress markers
   may indicate a root sampling artefact (wounding during harvest).

### Circadian considerations

Many plant genes (30-50% in Arabidopsis) show circadian oscillation.
Before clustering:
1. Run a periodogram (e.g., `MetaCycle` R package: JTK_CYCLE, Lomb-Scargle).
2. If a gene has significant periodicity (p < 0.05), consider analyzing
   circadian genes separately from non-circadian genes.
3. In rice/maize field experiments: diurnal temperature fluctuation
   can induce apparent circadian patterns that are actually
   temperature-driven. Flag these.

## Choosing Your Methods

### Smoothing: Spline vs LOESS

**Spline smoothing** (>= 5 time points):
- Fits smooth curves with continuous first and second derivatives
- Preserves inflection points (peak timing, transition points)
- df (degrees of freedom) selection: use generalized cross-validation (GCV)
- Risk: overfitting at the boundaries of the time course

**LOESS smoothing** (< 5 time points):
- Local polynomial regression, more robust with sparse data
- Span parameter controls smoothness (default 0.75; lower = more wiggly)
- Does not extrapolate beyond data range
- Risk: boundary bias at the first and last time points

### Clustering: Soft vs Hard

**Mfuzz (soft clustering, >= 1000 genes)**:
- Each gene has a membership score (0-1) for every cluster
- Captures genes with transitional or mixed patterns
- Useful when many genes show graded rather than categorical responses
- Cluster number selection: use `cselection()` or minimum centroid distance
- Plant-specific advantage: gene families often show subtly different
  paralog expression, captured by soft membership

**k-means (hard clustering, < 1000 genes)**:
- Each gene assigned to exactly one cluster
- Simple, fast, interpretable
- Cluster number selection: elbow method or gap statistic
- Works well for focused gene sets (e.g., DEG subset, TF family)

### Trend Analysis: STEM

- For >= 6 time points: fits predefined model profiles
- Each profile tests whether the number of genes assigned to it
  exceeds the expected number by chance
- Template profiles: monotonic up/down, transient up/down, biphasic
- Significance: permutation-based, corrected for multiple testing
- Output: colored significance plot showing enriched profiles

### Stage Specificity: Tau Index

- tau = 0: gene expressed uniformly across all tissues/stages
- tau = 1: gene expressed exclusively in one tissue/stage
- Calculated as: tau = sum(1 - xi/max(x)) / (n-1) where xi is expression in tissue i
- Plant reference thresholds:
  - tau > 0.8: tissue-specific
  - tau 0.5-0.8: tissue-enriched
  - tau < 0.5: broadly expressed
- Caveat: tau is sensitive to the set of tissues included. A gene that
  appears leaf-specific in a leaf-vs-root comparison may be broadly
  expressed when flower tissue is added.

## Common Pitfalls

### Batch effects masquerading as temporal patterns

If samples from the same time point were processed on different days
or by different people, batch effects can create spurious patterns.
Always check: is the pattern still present after batch correction?

### Normalization artifacts

DESeq2 normalization assumes most genes do not change. In dense time
courses with global expression shifts (e.g., senescence, heat shock),
this assumption may be violated. Check the distribution of normalized
counts per time point -- if the median shifts systematically, consider
using spike-in normalization or RUVseq.

### Over-interpreting cluster membership

A gene in cluster 3 is not necessarily "functionally related" to every
other gene in cluster 3. Always validate with:
1. GO enrichment per cluster
2. Known marker gene concordance
3. Promoter motif enrichment (plant-specific: check for known cis-elements
   in cluster members, e.g., G-box for light response, ABRE for ABA)

### Developmental asynchrony

In whole-organism or whole-organ sampling, developmental stages are
not perfectly synchronized across individuals. The expression pattern
you observe is a population average. A "gradual" increase may actually
be a sharp switch in individual plants that appears gradual due to
developmental asynchrony.

### Ignoring genotype effects

If your time course includes multiple genotypes (e.g., wild-type
plus mutant), cluster separately first, then compare. Pooling
genotypes can create artefactual clusters driven by genotype
differences rather than temporal dynamics.
