# Plant Phenotype Analysis Specifics

## Trial Design Types

| Design | When Used | Analysis Implication |
|--------|----------|---------------------|
| RCBD (Randomized Complete Block) | Standard field trial | Block as random effect |
| Alpha-Lattice / Lattice | Large trials with many entries | Row + Col as random |
| Augmented Design | Limited seed/space, unreplicated checks | Use checks for spatial correction |
| p-Rep Design | Partially replicated | Specialized analysis needed |

## Multi-Environment Trial Networks

- Connectedness: ensure genotypes are shared across environments for BLUP estimation
- Target Population of Environments (TPE): define which environments your variety targets
- Genotype × Environment (G×E): test for G×E significance before pooling environments

## Heritability Types

- **Broad-sense (H²)**: All genetic variance / total variance. Used for clonally propagated or inbred lines.
- **Narrow-sense (h²)**: Additive variance / total variance. Used for sexual breeding programs.
- **Plot-basis vs Entry-mean**: Entry-mean heritability is always higher — report which one you're using.
