# Multi-Environment Trial Summary

**Goal:** Summarize genotype performance across environments, assess stability
**Best for:** Multi-environment trials (≥2 environments)

## AMMI (Additive Main Effects and Multiplicative Interaction)

```r
library(agricolae)
ammi <- AMMI(env, genotype, rep, trait, data = pheno)
plot(ammi)
```

## GGE Biplot

```r
library(GGEBiplots)
gge <- GGEModel(pheno_matrix)
plot(gge, type = 1)  # Which-won-where pattern
plot(gge, type = 2)  # Mean vs stability
```

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| Number of IPC (AMMI) | 2 (IPC1, IPC2) | If IPC3 explains >10% of GxE variance, include it | Most GxE signal is captured in the first two components; additional PCs below 10% variance add noise without biological interpretation |
| GGE biplot type | type=1 (which-won-where) | Use type=2 for mean-vs-stability; type=3 for ranking genotypes; type=4 for discriminating power vs representativeness | Different biplot types answer fundamentally different breeding questions — choose based on your goal, not the default |
| Minimum environments | 3 | For AMMI with reliable IPC estimation, >=6 environments strongly preferred | With <6 environments, the residual GxE matrix has few degrees of freedom and IPC estimates become unstable |
| Scaling (GGE) | environment-centered (column metric preserving) | Use genotype-centered if the focus is on genotype relationships rather than GxE | Scaling determines whether the biplot emphasizes GxE interaction (env-centered) or genotype similarity (genotype-centered) |
| Missing value handling | Impute via EM-AMMI | If >10% missing data, remove genotypes/environments with excessive missingness | EM-AMMI imputation is reliable for sparse missing values but breaks down with systematic missing data patterns |

## Plant-Specific Notes
- AMMI/GGE help identify broadly adapted vs specifically adapted genotypes
- "Which-won-where" biplot shows which genotype performs best in each environment
- For breeding: select genotypes stable across target environments
- **MET with >3 environments**: AMMI requires at least 3 environments to estimate the first interaction PC. However, for publication-quality GxE analysis, aim for >=6 environments with >=2 replicates each. Fewer environments limit the detectable GxE patterns and produce unreliable mega-environment classifications.
- **Cross-year interpretation**: Environment is a combination of location + year. Two trials at the same location in different years are different environments. If you have multi-year data, include year as an environmental factor — GxE patterns can shift dramatically between years due to rainfall and temperature differences.
- **Heritability check**: Before interpreting AMMI/GGE results, check the within-environment heritability (H2). Environments with H2 < 0.3 have unreliable genotype rankings and will distort GGE biplots. Consider excluding low-heritability environments from the analysis.

## Common Errors
| Error | Cause | Solution |
|-------|-------|----------|
| "Cannot compute AMMI" | Missing values in phenotype matrix | Impute or remove genotypes with missing data |
| AMMI IPC1 explains >80% of GxE | GxE is dominated by a single pattern (e.g., one environment very different from others) | Report this as a finding; the dominant IPC may represent a single outlier environment — check if that environment should be analyzed separately |
| GGE biplot "which won where" shows no clear winners | GxE is small relative to main genotype effect | Use mean-vs-stability plot (type=2) instead; which-won-where is uninformative when GxE is weak |
| Genotype rankings differ between AMMI and GGE | AMMI includes main additive effects in the model; GGE merges G+GxE | Report both; AMMI ranks are more relevant for broad adaptation; GGE ranks are better for specific adaptation |
| Negative IPC variance estimates | Too few environments or too many missing values | Reduce the number of IPCs; check for environments with near-zero variance |
