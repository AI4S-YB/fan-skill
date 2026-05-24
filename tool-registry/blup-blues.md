# BLUP / BLUE Estimation — Genotype Effect Prediction for Plant Breeding

**Goal:** Estimate Best Linear Unbiased Predictions (BLUP) and Best Linear Unbiased Estimates (BLUE) for genotype effects from field trial data
**Approach:** Mixed linear models: genotype as random (BLUP, shrinks effects toward mean) or fixed (BLUE, direct estimates) → extract predicted genotype values → optionally estimate heritability (H^2)
**Best for:** Ranking genotypes for selection, preparing genomic selection training data, multi-environment trial (MET) analysis, and calculating genetic gain

## Prerequisites
- R 4.0+, lme4 (for simple models) or sommer (for complex variance structures)
- Phenotype data with genotype identifiers, environment identifiers (for MET), and optional block/replicate information
- Balanced or unbalanced designs — mixed models handle both, but convergence is easier with balanced data
- For heritability: replicate observations per genotype (at minimum, within-environment replication)

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(sommer)  # Preferred for plant breeding: supports complex variance structures

# ── Option A: Single-Environment BLUP (genotype as random) ──
# Use when ranking genotypes within one trial or extracting GS training values
ans_blup <- mmer(
  fixed  = ${TRAIT} ~ ${FIXED_EFFECTS},    # e.g., ~ 1 (intercept only) or ~ Rep
  random = ~ vs(${GENOTYPE_COL}, Gu = ${K_MATRIX}),  # Gu = kinship (optional, for GBLUP)
  rcov   = ~ units,
  data   = ${PHENO_DATA}
)
blups <- ans_blup$U$`${GENOTYPE_COL}`$${TRAIT}
# Reliability (r²): 1 - (PEV / genetic variance), where PEV from summary(ans_blup)

# ── Option B: Single-Environment BLUE (genotype as fixed) ──
# Use when you need unbiased estimates (e.g., comparing specific genotypes, not ranking)
ans_blue <- mmer(
  fixed  = ${TRAIT} ~ ${GENOTYPE_COL} + ${FIXED_EFFECTS},
  random = ~ ${RANDOM_EFFECTS},            # e.g., ~ block if RCBD
  rcov   = ~ units,
  data   = ${PHENO_DATA}
)
blues <- summary(blue_model)$beta

# ── Option C: Multi-Environment BLUP (genotype × environment) ──
# Use for MET analysis: extracts stable genotype effects across environments
ans_met <- mmer(
  fixed  = ${TRAIT} ~ ${ENV_COL},
  random = ~ vs(${GENOTYPE_COL}) + vs(${ENV_COL}:${GENOTYPE_COL}),
  rcov   = ~ vs(${ENV_COL}),              # Heterogeneous error variance per environment
  data   = ${PHENO_DATA}
)
blups_main      <- ans_met$U$`${GENOTYPE_COL}`$${TRAIT}       # Main genotype effect
blups_interact  <- ans_met$U$`:${ENV_COL}:${GENOTYPE_COL}`$${TRAIT}  # G×E interaction

# ── Option D: Heritability Estimation ──
# Broad-sense heritability: H² = Vg / (Vg + Vge/n_env + Ve/(n_env * n_rep))
# Using sommer variance components:
vc <- summary(ans_met)$varcomp
Vg  <- vc[grep("${GENOTYPE_COL}", rownames(vc)), "VarComp"]
Vge <- vc[grep("${ENV_COL}:${GENOTYPE_COL}", rownames(vc)), "VarComp"]
Ve  <- vc[grep("units", rownames(vc)), "VarComp"]
n_env <- length(unique(${PHENO_DATA}$${ENV_COL}))
n_rep <- ${N_REPS}  # Harmonic mean for unbalanced designs
H2   <- Vg / (Vg + Vge/n_env + Ve/(n_env * n_rep))

# ── Alternative: lme4 for simpler models (faster, less flexible) ──
library(lme4)
blup_lme4 <- lmer(${TRAIT} ~ ${FIXED} + (1|${GENOTYPE_COL}) + (1|${ENV_COL}:${GENOTYPE_COL}),
                  data = ${PHENO_DATA})
blups <- ranef(blup_lme4)$${GENOTYPE_COL}
```

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| Genotype as random (BLUP) | Yes | When estimating fixed genotypic values (BLUE) for specific comparisons → Fixed | BLUP shrinks toward mean — better for ranking, worse for absolute values |
| Genotype as fixed (BLUE) | No | When you have few genotypes (<10) and want unbiased estimates → Fixed | With few levels, random effect variance is poorly estimated |
| Include G×E in random | Yes (MET) | Single-environment trial → Not applicable | G×E is confounded with residual in single-env models |
| `Gu = kinship` | NULL | When using genomic relationships (GBLUP) → provide additive relationship matrix | Improves accuracy when pedigree or markers are available |
| Heterogeneous error (`rcov = ~ vs(env)`) | Yes (MET) | When environments have very similar error variance → `~ units` | Heterogeneous variance is more realistic, but simpler model converges faster |
| Spatial correction | Not included | When trial shows spatial trends (row/column effects) → add `~ vs(row) + vs(col)` or spline terms | Ignoring spatial variation inflates error variance and reduces H² |
| `n_rep` for H² | observed reps | Unbalanced MET → use harmonic mean of reps | Arithmetic mean overestimates H² in unbalanced designs |

---

## Plant-Specific Notes

### MET Design Types and Model Implications

| Design | Typical crops | Model recommendation |
|--------|--------------|---------------------|
| **RCBD (Randomized Complete Block)** | Most field crops (maize, wheat, rice, soybean) | Block as fixed or random within environment; genotype as random |
| **Alpha-Lattice / Incomplete Block** | Large breeding trials (>=100 genotypes) | Incomplete block nested within rep as random |
| **Augmented Design** | Preliminary yield trials, GRIN/GenBank evaluations | Checks (repeated varieties) as fixed; test genotypes as random |
| **Row-Column Design** | Forestry, horticulture, spaced-plant trials | Row + Column as random; add spatial AR1 x AR1 correlation if available |
| **p-rep Design (partially replicated)** | Early-generation selection (F3-F5) | Unreplicated genotypes borrow strength from replicated checks; use spatial correction |
| **Perennial crop repeated measures** | Apple, grape, tea, oil palm | Add year and genotype×year as random; account for autocorrelation across years |

### Spatial Correction Integration
- Row and column effects should be modeled as random (not fixed) in plant breeding trials
- For large trials (>1000 plots), 2D P-spline spatial correction (`SpATS` package in R) significantly improves H²
- Spatial trend often accounts for 20-40% of residual variance — check variogram or heatmap of residuals before deciding
- **When to add spatial terms**: residual variogram shows spatial autocorrelation or field heatmap reveals fertility gradients

### Polyploid species
- For autopolyploids (potato, alfalfa, sugarcane): additive relationship matrix must account for double-reduction and allele dosage
- For allopolyploids (wheat, cotton, canola): use subgenome-specific kinship matrices as additional random effects to partition subgenome contributions
- Genomic BLUPs (GBLUP) from markers may outperform pedigree-based BLUPs for polyploids where pedigree records are complex

### Cross-pollinated vs self-pollinated crops
- **Self-pollinated (rice, wheat, soybean, tomato)**: genotypes are nearly homozygous lines — BLUPs are stable across generations
- **Cross-pollinated (maize, sunflower, rye)**: hybrid or open-pollinated varieties — use testcross BLUP (genotype effect conditional on tester); do NOT use per se BLUP for hybrid prediction

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Model failed to converge" | Unbalanced design, too many random effects, or insufficient replication | Simplify random structure; remove higher-order interactions; try `lme4` with different optimizers (`Nelder_Mead`, `bobyqa`) |
| "Singular fit" in lme4 | Random effect variance estimated as zero (e.g., G×E = 0) | Remove the zero-variance term from the model; report that the effect was not detectable |
| BLUPs all near zero with large standard errors | Too few replicates per genotype (<2) | BLUP shrinkage is severe with little data — increase replication or use spatial correction to borrow information |
| H² estimate >1.0 or negative | Wrong variance component partitioning; outlier genotypes inflating Vg | Check for outliers (studentized residuals > 3); verify model formula; negative H² → set to 0 |
| BLUPs change dramatically with a few added genotypes | Unbalanced data with influential genotypes | Use robust BLUP or check leverage; ensure genotype connectivity across environments |
| G×E variance component zero | Genotypes rank consistently across environments (no crossover interaction) | Not an error — report it; this is good news for breeding (broad adaptation) |
| H² very different between BLUP and ANOVA-based method | BLUP accounts for unbalanced data and shrinkage; ANOVA does not | Always report BLUP-based H² for unbalanced MET data |
| "contrasts can be applied only to factors with 2 or more levels" | A factor has only one level after subsetting | Check that you have >=2 environments, >=2 genotypes, and replication per level |
| Lost degrees of freedom for fixed genotype model | Too many genotypes as fixed effects + small residual df | Switch to random genotype (BLUP); with >50 genotypes, fixed model is wasteful |
