# Heritability Estimation — Broad-Sense & Narrow-Sense

**Goal:** Estimate heritability from replicated phenotype data to quantify genetic vs environmental variance — essential before GWAS or GS
**Approach:** Fit mixed linear model with genotype as random effect → partition variance components → H² = Vg / (Vg + Ve/r)
**Best for:** All quantitative traits; mandatory quality-control step before investing in genomic analysis

## Prerequisites
- R 4.0+, lme4 (REML), sommer (flexible variance structures), lmerTest (p-values)
- Phenotype data: genotype/line IDs, replication structure, environment labels
- For narrow-sense heritability: pedigree or genomic relationship matrix

---

## Code Skeleton

This is NOT a fixed script. Adapt `${PLACEHOLDERS}` to your data.

```r
library(lme4)
library(lmerTest)       # optional: for p-values on fixed effects
# Alternative: library(sommer)  # uncomment when needed (see When lme4 vs sommer)

# ── Step 1: Load phenotype data ──
pheno <- ${PHENO_DATA}  # data.frame with columns: genotype, trait, environment, rep, block
str(pheno)

# ── Step 2: Single-environment broad-sense heritability (H²) ──
# Model: trait = μ + genotype(random) + rep(fixed) + ε
model_se <- lmer(
  ${TRAIT} ~ ${REP} + (1 | ${GENOTYPE}),
  data = pheno,
  REML = TRUE
)
vc <- as.data.frame(VarCorr(model_se))
Vg <- vc$vcov[1]                     # genetic variance (first random term)
Ve <- attr(VarCorr(model_se), "sc")^2  # residual variance
n_reps <- ${N_REPS}                   # number of replicates

H2 <- Vg / (Vg + Ve / n_reps)        # broad-sense heritability (line-mean basis)
cat(sprintf("Single-environment H² = %.3f\n", H2))

# ── Step 3: Multi-environment broad-sense heritability (MET H²) ──
# Model: trait = μ + genotype + environment + genotype×env + rep(env) + ε
# All random except intercept — this is the standard MET model
model_met <- lmer(
  ${TRAIT} ~ (1 | ${GENOTYPE}) + (1 | ${ENVIRONMENT}) +
             (1 | ${GENOTYPE}:${ENVIRONMENT}) + (1 | ${ENVIRONMENT}:${REP}),
  data = pheno,
  control = lmerControl(optimizer = "bobyqa")
)

vc_met <- as.data.frame(VarCorr(model_met))
Vg_met  <- vc_met$vcov[vc_met$grp == "${GENOTYPE}"]
Vge     <- vc_met$vcov[vc_met$grp == "${GENOTYPE}:${ENVIRONMENT}"]
Ve_met  <- attr(VarCorr(model_met), "sc")^2
n_env   <- ${N_ENVIRONMENTS}
n_rep   <- ${N_REPS}

# Line-mean heritability across environments
H2_met <- Vg_met / (Vg_met + Vge / n_env + Ve_met / (n_env * n_rep))
cat(sprintf("Multi-environment H² = %.3f\n", H2_met))

# ── Step 4: Narrow-sense heritability (h²) via sommer ──
# Requires pedigree or genomic relationship matrix (A or G)
# Uncomment when relationship matrix available:
# library(sommer)
# model_ns <- mmer(
#   fixed = ${TRAIT} ~ 1,
#   random = ~ vs(${GENOTYPE}, Gu = ${A_MATRIX}),
#   data = pheno
# )
# Va <- summary(model_ns)$varcomp[1, "VarComp"]    # additive genetic variance
# Ve_ns <- summary(model_ns)$varcomp[2, "VarComp"]  # residual variance
# h2 <- Va / (Va + Ve_ns / n_reps)
# cat(sprintf("Narrow-sense h² = %.3f\n", h2))

# ── Step 5: Report variance components ──
# Always report ALL components, not just H²
cat("\n─── Variance Components Report ───\n")
cat(sprintf("Vg (genetic):         %.4f\n", Vg))
cat(sprintf("Ve (residual):         %.4f\n", Ve))
cat(sprintf("H² (broad-sense):      %.3f\n", H2))
# With MET: also report Vge (G×E interaction variance)
```

---

## Broad-Sense (H²) vs Narrow-Sense (h²)

| Measure | Definition | What it captures | Used for |
|--------|------|------|----------|
| **H²** (broad-sense) | Vg / Vp | All genetic variance (additive + dominance + epistatic) | Clonal propagation, inbred lines, trait screening |
| **h²** (narrow-sense) | Va / Vp | Additive genetic variance only | Genomic selection, breeding value prediction, selection response |
| Relationship | h² <= H² | — | Always ≤ H²; equal only when dominance = 0 (inbreds at fixation) |

---

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| `REML` | TRUE | Fixed effects comparison → use ML temporarily | REML gives unbiased variance estimates |
| `optimizer` | bobyqa | Convergence warning → try "Nelder_Mead" or "nloptwrap" | Different optimizers handle boundary constraints differently |
| `Model type` | `lmer` (lme4) | Spatial correction, heterogeneous variance needed → `mmer` (sommer) | sommer supports flexible covariance structures |
| `Genotype term` | `(1\|genotype)` random | Genotypes are fixed varieties (checks) → fixed effect | Inference goal determines random vs fixed |
| `Replicate term` | `(1\|rep)` or `(1\|env:rep)` | Rep is unreplicated within env → omit | To avoid non-identifiable parameter |

---

## When to Use lme4 vs sommer

| Scenario | Package | Rationale |
|----------|:---:|------|
| Simple H² estimation, balanced design | **lme4** | Simpler syntax, faster, well-tested |
| Multi-environment with heterogeneous error variance | **sommer** | `mmer()` allows `rcov = ~ vs(units)` for per-environment residual variance |
| Spatial field trial (row × column effects) | **sommer** | `spl2Da()` and `spl2Db()` for 2D spatial splines |
| Pedigree-based h² (with A matrix) | **sommer** | `vs(genotype, Gu = A)` native; lme4 requires manual Cholesky decomposition |
| Genomic h² (with G matrix) | **sommer** or **BGLR** | sommer: `vs(genotype, Gu = G)`; BGLR: RKHS model |
| Non-Gaussian traits (binary, count) | **MCMCglmm** or **brms** | lme4 `glmer` is an option but Bayesian methods handle boundary issues better |

---

## Plant-Specific Notes

### Multi-Environment Trials (MET) — the plant-breeding norm
- **Always use MET when available.** Single-environment H² is environment-specific and NOT generalizable to other environments or years.
- The MET model `(1|genotype) + (1|env) + (1|genotype:env) + (1|env:rep)` is the standard for plant breeding programs.
- G×E variance (Vge) is typically >0 for yield and stress-tolerance traits — if Vge = 0, check the data; it may be a convergence artifact with few environments.
- **Minimum environments:** 3+ environments for reliable H². With 2 environments, G×E is estimated but unreliable.

### Spatial Field Trial Correction — essential for large field plots
- Row and column effects from soil gradients can inflate Ve, deflating H².
- Use sommer with `spl2Da(row, col)` or `spl2Db(row, col)` for 2D spatial spline corrections.
- Block effects capture only coarse spatial variation; splines capture fine-scale gradients.
- Example check: plot residuals by field position. Diagonal striping or edge effects signal uncorrected spatial variation.

### Self-pollinated crops (rice, wheat, soybean)
- Inbred lines are homozygous: H² (broad-sense) captures fixed line differences.
- For elite lines from same crossing program, H² may be low simply because genetic variance is low — not a model failure.
- Narrow-sense h² within a biparental population derived from two inbred parents is well-defined (RILs, DH lines).

### Outcrossing / clonal species
- Clonally propagated crops (potato, sugarcane, fruit trees): H² is directly relevant since clones capture total genetic value.
- For hybrid crops (maize): h² within heterotic groups matters for breeding decisions; H² across groups is less interpretable.
- Perennial fruit trees: spatial effects from orchard topography are critical — always run the spatial model.

---

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| H² > 1 | Negative variance component estimate (boundary constraint) | Use REML; check design balance; consider Bayesian model (MCMCglmm) |
| H² ≈ 0 or negative | Low genetic variance, high residual error, or design confounded | May be real — calculate CV%; if design is balanced, report honestly |
| "boundary (singular) fit" | Variance component estimated at zero | Simplify random structure; check if genotype term has enough levels |
| "Model failed to converge" with lmer | Optimizer stuck, bad scaling | Scale trait: `scale(trait)`; change optimizer to "Nelder_Mead" or "bobyqa" |
| G×E variance estimated as zero | Too few environments (<3); or environments are very similar | Minimum 3 distinct environments; report uncertainty |
| H² differs wildly between environments | Genuine G×E — heritability IS environment-dependent | Report per-environment H² + overall MET H²; never average per-env H² |
| Standard error of H² is huge | Too few genotypes (<30) or unbalanced data | Increase genotype count; report SE alongside H² estimate |
| Negative Vg from sommer | Boundary hit during REML; parameter at edge of feasible space | Constrain variance with `bounds` argument; or use `EMMA` algorithm |
