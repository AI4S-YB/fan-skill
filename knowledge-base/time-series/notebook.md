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

### Sampling frequency: the Nyquist of plant biology

Choosing how often to sample is a trade-off between biological resolution
and practical constraints (cost, labor, tissue availability).

**Sampling frequency by time scale**:

| Biological process | Recommended interval | Minimum to capture peak |
|-------------------|---------------------|------------------------|
| Circadian rhythm | 2-4 hr | 4 hr (Nyquist for 24 hr period) |
| Rapid stress (ROS, Ca²⁺) | 5-15 min | 15 min |
| Transcriptional stress response | 0.5-2 hr | 2 hr |
| Hormone response | 1-4 hr | 4 hr |
| Developmental (days) | 1-3 days | 3 days |
| Seed maturation | 2-5 days | 5 days (process takes weeks) |

**The "missed peak" problem**:
- If your sampling interval is 4 hr and the true expression peak lasts 2 hr,
  you have a 50% chance of missing it entirely.
- Solution for rapid responses: pilot experiment with denser sampling on
  a few marker genes (qRT-PCR) to determine the peak window, then design
  the full RNA-seq time course around that window.
- Multi-phasic responses (e.g., early peak at 1 hr, late peak at 24 hr)
  require log-spaced time points: 0, 0.5, 1, 2, 4, 8, 24 hr covers both
  phases efficiently.

**Continuous vs discrete sampling**:
- Destructive sampling (most plant RNA-seq): you harvest a whole leaf/root,
  so each time point is an independent biological sample. This inflates
  biological variability. Compensate with higher replicate numbers (n >= 3
  per time point).
- Non-destructive sampling is rare in plants but possible for certain systems
  (e.g., leaf punch on large leaves, or imaging-based reporters). If available,
  it dramatically reduces inter-sample variability.

### Tissue atlas interpretation

Plant tissue atlases (e.g., rice eFP browser, Arabidopsis TraVA, maize
Expression Atlas) provide spatial context. When interpreting clusters:

1. Cross-reference cluster members with known tissue markers.
2. A cluster enriched in "photosynthesis" GO terms that peaks at midday
   in leaf tissue is expected -- do not overinterpret.
3. Clusters where root-specific genes co-express with stress markers
   may indicate a root sampling artefact (wounding during harvest).

### Plant developmental gradients: more than just time

Plant development is not purely chronological -- it has a spatial axis
that is uniquely important compared to animal systems:

**Leaf age gradients**:
- In monocots (rice, maize, wheat): leaves emerge sequentially, with leaf 1
  being the oldest and the flag leaf the youngest. Sampling "leaf tissue" at
  different positions along the shoot is sampling a developmental gradient,
  not a single tissue type.
- The leaf age gradient is a continuum: young (sink) leaves import carbon,
  mature (source) leaves export carbon. Expression profiles reflect this
  metabolic transition.
- Recommendation: for developmental studies in grasses, record leaf number
  AND leaf position. A metadata field "leaf_position" (flag leaf = 1,
  flag-1 = 2, etc.) enables post-hoc parsing of true developmental signal
  from leaf-age confounding.

**Fruit ripening gradients**:
- In fleshy fruits (tomato, grape, citrus): ripening proceeds through
  color stages (mature green → breaker → turning → red ripe). These are
  discrete, visually scorable stages that correlate with massive
  transcriptomic and metabolomic shifts.
- Important: ripening is often asynchronous within a cluster (truss).
  Sampling "breaker-stage" fruit from different trusses without controlling
  for truss position introduces variability. Specify truss number in metadata.
- Ethylene climacteric (tomato, apple, banana) vs non-climacteric (grape,
  citrus, strawberry): the transcriptional cascade differs fundamentally.
  Do not pool these species types in comparative analyses.

**Seed maturation gradients**:
- Developing seeds progress through: embryogenesis → grain filling →
  desiccation → dormancy acquisition. Each phase has distinct transcriptional
  programs.
- In cereals (rice, wheat, barley): days after anthesis (DAA) is the
  standard staging system. Endosperm and embryo should be profiled
  separately -- their transcriptomes diverge dramatically after 5-10 DAA.
- In legumes (soybean, pea): seed size rather than DAA is often the more
  reliable staging metric because flowering is indeterminate (new flowers
  appear while pods are already filling), so DAA mixes different
  developmental stages on the same plant.

### Circadian considerations

Many plant genes (30-50% in Arabidopsis) show circadian oscillation.
Before clustering:
1. Run a periodogram (e.g., `MetaCycle` R package: JTK_CYCLE, Lomb-Scargle).
2. If a gene has significant periodicity (p < 0.05), consider analyzing
   circadian genes separately from non-circadian genes.
3. In rice/maize field experiments: diurnal temperature fluctuation
   can induce apparent circadian patterns that are actually
   temperature-driven. Flag these.

### Circadian experimental design: avoiding systematic confounding

Circadian rhythms are the most pervasive source of temporal confounding
in plant expression studies. Poor design can render your entire experiment
uninterpretable.

**The "time-of-day" trap**:
- If you sample "control" plants at 9 AM and "treated" plants at 3 PM, any
  differential expression could be treatment effect OR circadian effect.
  You cannot distinguish them.
- Solution: always collect control and treated samples at matched Zeitgeber
  time (ZT). If treatment application itself takes time, randomize the order
  of control vs treated harvesting within each time point.

**Free-running vs entrained conditions**:
- Entrained (light-dark cycles): appropriate for most experiments. Plants
  are in their natural rhythm. But the light-dark transition itself causes
  transcriptional bursts that are not "circadian" (they are light-responsive).
- Free-running (constant light or constant dark): isolates the endogenous
  circadian clock from light-driven responses. Essential if your biological
  question is specifically about the circadian clock mechanism.
- Important: after 48+ hr in constant light, plants show metabolic stress
  (sugar accumulation, photoinhibition). Limit free-running experiments to
  48-72 hr unless the clock mechanism itself is the target.

**Temperature cycles as a confounder**:
- In field experiments, temperature varies with time of day. Many
  "circadian" genes are actually temperature-responsive (thermocycle).
  True circadian genes continue oscillating under constant temperature;
  thermocycle-driven genes do not.
- In growth chamber experiments, verify that your temperature control is
  actually stable. Cheap growth chambers can drift 2-3°C diurnally,
  enough to induce thermocycle-responsive genes.

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

### Soft vs Hard Clustering: Decision Framework

Choosing between Mfuzz and k-means is not just about gene count.
Consider the underlying biology:

**Choose Mfuzz when**:
1. Your time points represent a continuous developmental gradient (not
   discrete stages). Leaf aging, fruit ripening, and stress progression
   are continuous processes -- genes transition gradually between states,
   and soft membership captures this continuum.
2. You have paralog-rich gene families. In polyploid plants (wheat, cotton,
   oilseed rape), homeologs often show subtly different expression timing.
   Hard clustering forces them apart; Mfuzz assigns them to neighboring
   clusters with overlapping membership.
3. The biological question is about regulatory gradients rather than
   discrete states. Mfuzz's membership values can be interpreted as
   "strength of participation" in a temporal program.
4. Your data is inherently noisy (field-collected samples, variable
   developmental staging). Mfuzz's soft assignment is more robust to noise
   because it doesn't force ambiguous genes into a single cluster.

**Choose k-means when**:
1. Your time points represent discrete, well-separated conditions (e.g.,
   four different tissues, not really a time course, or clear before/after
   treatment). Hard assignment is appropriate.
2. You need simple, interpretable output for a non-computational audience.
   "Gene X belongs to cluster 3" is easier to communicate than "Gene X has
   membership 0.6 in cluster 3 and 0.3 in cluster 4."
3. You are working with a focused candidate gene set (< 500 genes). Soft
   clustering with too few genes produces unstable membership values.
4. You plan to do downstream functional validation on cluster representatives.
   Hard clusters give you clean gene lists for GO enrichment.

**Hybrid approach (recommended for large studies)**:
1. Run both Mfuzz and k-means on the same data
2. For each gene, compare its hard cluster assignment vs its primary
   (highest-membership) Mfuzz cluster
3. Genes with concordant assignment (> 80% typically) are "core" cluster
   members -- high confidence. Focus functional interpretation on these.
4. Genes with discordant assignment are "transitional" -- interesting
   biology (e.g., bifunctional regulators) but lower confidence. Report
   separately.

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

### Stage-Specificity Metrics: Beyond Tau

The tau index is the most widely used metric, but it has limitations.
Consider complementary metrics:

**Tau limitations in plant contexts**:
- Tau assumes all tissues/stages are equally "different" from each other.
  In a developmental series (leaf 1→leaf 8), leaf 3 and leaf 4 are more
  similar than leaf 1 and leaf 8. Tau treats them as equally distinct.
- Tau is sensitive to outliers. A single time point with abnormally high
  expression (e.g., a stress-responsive spike in one sample) inflates tau
  even if the gene is broadly expressed at all other stages.
- For circadian genes, tau computed across a diurnal time course is
  misleading -- a gene peaking at ZT4 and ZT28 looks "specific" to two
  time points, but it's actually rhythmic.

**Complementary metrics**:
- **Tau with temporal weighting**: modify the tau formula to weight
  neighboring time points more similarly. For developmental gradients,
  use `tau_w = sum(w_i * (1 - xi/max(x))) / (n-1)` where w_i is a
  distance-based weight.
- **Specificity score (SP)**: SP = (expression in tissue i) / (sum of
  expression across all tissues). Unlike tau, SP is per-tissue, not
  per-gene. Useful for identifying tissue-specific genes when you have
  clear tissue labels.
- **JS distance (Jensen-Shannon)**: measures how different a gene's
  expression distribution is from a uniform distribution. Less sensitive
  to outliers than tau.
- **Practical recommendation**: report tau + one complementary metric
  (SP or JS). If both agree, the specificity call is robust. If they
  disagree, investigate the gene's expression pattern manually.

### Stress Response Time-Course Design

Plant stress experiments have unique design challenges that differ from
developmental or circadian studies:

**Acute vs chronic stress**:
- Acute stress (sudden drought, heat shock, wounding): triggers a rapid
  transcriptional cascade (minutes to hours). Early time points (0-4 hr)
  are the most informative window. Dense sampling early: 0, 0.5, 1, 2, 4 hr.
- Chronic stress (gradual drought, progressive nutrient deficiency, prolonged
  cold): involves acclimation and metabolic reprogramming. Sampling should
  extend to days or weeks. Log-spaced time points work best.
- Many genes show a "spike then decline" pattern in acute stress but a
  "gradual ramp" in chronic stress. If you only sample one time regime,
  you mischaracterize the gene's response.

**Recovery time points**:
- Adding a "recovery" time point (re-water after drought, return to normal
  temperature after heat) reveals which stress responses are reversible
  and which represent lasting reprogramming.
- Recovery-responsive genes (return to baseline within hours) vs
  memory-responsive genes (maintain altered expression for days) have
  different biological significance. Recovery genes are involved in acute
  homeostasis; memory genes in epigenetic stress priming.
- Recommend at least 2 recovery time points: one early (1-2 hr post-recovery)
  and one late (24-48 hr post-recovery).

**Replicate design in stress experiments**:
- Stress responses are inherently more variable than developmental
  trajectories. Individual plants may differ in how quickly they perceive
  and respond to stress.
- Minimum 3 biological replicates per time point for stress studies (vs
  2 for tightly controlled developmental experiments).
- If possible, use paired design: sample the same plant before and after
  stress. Not always feasible (destructive sampling of roots/leaves), but
  when available (e.g., leaf punch on large-leaf species), it dramatically
  reduces inter-plant variability.

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

### Circadian confounding in non-circadian experiments

This is one of the most insidious pitfalls in plant time-series analysis.
Even if you are not studying circadian rhythms, the circadian clock is
confounding your data:

**How to detect circadian confounding**:
1. If your experiment lasted > 6 hours, circadian effects are present.
2. Check known clock genes (CCA1, LHY, TOC1, PRR family, GI in
   Arabidopsis; OsCCA1, OsPRR1, OsGI in rice; ZmCCA1, ZmTOC1 in maize).
   If these show rhythmic patterns in your data, circadian effects are
   real and pervasive.
3. A cluster with ~24 hr oscillation that appears in "control" samples
   but disappears in "treated" samples may indicate treatment-induced
   clock disruption -- a biologically interesting finding, not a nuisance.

**How to mitigate**:
1. If possible, collect all samples at the same Zeitgeber time (ZT).
   This eliminates circadian variation but limits your experiment length
   to 24-hr multiples.
2. If you must sample across different ZTs, include ZT as a covariate
   in your statistical model: `expression ~ treatment + ZT + treatment:ZT`.
   The ZT term absorbs circadian variation.
3. In differential expression analysis between time-matched treatment
   and control samples, circadian effects cancel out (assuming the
   treatment does not alter the clock).
4. Never compare a treatment sample collected at ZT6 to a control sample
   collected at ZT12 and attribute the difference to treatment.

### Ignoring genotype effects

If your time course includes multiple genotypes (e.g., wild-type
plus mutant), cluster separately first, then compare. Pooling
genotypes can create artefactual clusters driven by genotype
differences rather than temporal dynamics.

### Developmental asynchrony: deeper treatment

Developmental asynchrony is pervasive in plant studies and can create
misleading "gradual" patterns:

**Quantifying asynchrony**:
- For a developmental stage defined by a morphological marker (e.g.,
  "anthesis" in rice), sample 20+ individuals and score the degree of
  variation. If anthesis spans a 3-day window in your population, your
  "0 DAA" sample actually represents 0-3 DAA.
- Impact: a gene with a sharp, 6-hour expression peak at true 0 DAA will
  appear to have a broad, 3-day peak in pooled samples.
- Mitigation: tag individual flowers/panicles and harvest based on
  individual developmental timing, not population average.

**When asynchrony is informative**:
- In some cases, developmental asynchrony is itself the phenotype.
  For example, drought escape in wheat: early-flowering genotypes
  "escape" terminal drought. The expression difference between
  early- and late-flowering genotypes at the same calendar date
  reflects both developmental stage and genotype effects.
- Do not "correct" for this asynchrony -- it's the biological signal.
  Instead, report both calendar-date and developmental-stage analyses.
