# Plant Enviromics Specific Considerations

## Why Plant Enviromics Differs from Standard G×E Analysis

| Aspect | Standard G×E Analysis | Plant Enviromics |
|--------|----------------------|------------------|
| Environmental data | Environment as categorical label | Quantitative environmental covariates (climate, soil) |
| G×E interpretation | Statistical significance | Biological mechanism (which variable drives interaction) |
| Scale | Single experiment | Multi-year, multi-location breeding trials |
| Goal | Detect G×E | Predict in new environments; delineate TPE |
| Temporal dimension | Ignored | Year-to-year climate variation explicitly modeled |
| Spatial dimension | Trial locations as points | Continuous environmental gradients |

## Defining Environmental Covariates

### What makes a good environmental covariate?

A useful environmental covariate should be:
1. **Biologically relevant**: Tied to crop physiology (not arbitrary)
2. **Measurable/derivable**: From weather stations, satellites, or soil maps
3. **Variable across environments**: A constant has no explanatory power
4. **Independent enough**: High collinearity reduces interpretability

### Common covariate categories for plants

| Category | Example Variables | Relevance |
|----------|------------------|-----------|
| Temperature | GDD, Tmean, Tmax > 30C count, Tmin < 0C count | Growth rate, heat/cold stress |
| Moisture | Total precipitation, rainy days, VPD, soil moisture | Water availability, drought stress |
| Radiation | Solar radiation, PAR, sunshine hours, photoperiod | Photosynthesis, photoperiod sensitivity |
| Soil | pH, organic matter, CEC, available P/K, texture | Nutrient availability, rooting depth |
| Management | Sowing date, density, N application, irrigation | Agronomic practices |

### Window-based covariates

The most biologically meaningful covariates are computed over specific
phenological windows, not the whole growing season:

| Crop | Critical Window | Duration | Key Variable |
|------|----------------|----------|--------------|
| Wheat | Anthesis +/- 10 days | 20 days | Tmax (heat stress at flowering) |
| Wheat | Grain filling | ~30 days | Tmean, VPD |
| Maize | Silking +/- 15 days | 30 days | Precip, soil moisture |
| Rice | Panicle initiation to heading | ~30 days | Tmin (cold-induced sterility) |
| Soybean | Pod filling (R5-R6) | ~20 days | Precip, soil moisture |
| Cotton | Boll development | ~40 days | Tmean, solar radiation |

### When to compute what?

- **If you have daily weather**: Compute window-specific GDD, stress counts,
  VPD, cumulative precipitation
- **If you only have monthly averages**: Use growing-season means;
  acknowledge the coarse resolution in interpretation
- **If you only have WorldClim bioclimatic variables**: Use BIO1 (annual
  mean temp), BIO12 (annual precip), BIO5 (max temp of warmest month),
  BIO15 (precip seasonality); note these are 30-year averages, not
  year-specific

## Multi-Environment Trial (MET) Design in Plant Breeding

### Common MET structures

1. **Incomplete block within location**: Standard RCBD / alpha-lattice
2. **Sparse testing**: Not all genotypes in all environments (genomic
   prediction can fill gaps)
3. **Phenotyping networks**: Regional variety trials (e.g., uniform
   nurseries)
4. **Managed stress trials**: Controlled drought, heat, or nutrient stress
   at the same location but different treatments

### Design implications for enviromics analysis

- **Sparse testing**: FA models can still work but genetic correlations
  between poorly connected environments have large SE
- **Managed stress**: The "treatment" is part of the environment definition.
  Treat managed stress environments as separate but related environments
- **On-farm trials**: Often high CV and unbalanced. Report per-environment
  heritability; low-heritability environments should be downweighted

## G×E in the Context of Climate Change

### Why enviromics matters more than ever

- Historical trial data embeds genotype responses to past climate variation
- Enviromics links this variation to climate variables, enabling
  prediction under future climates
- Key question: "Which current TPE is most similar to the future climate
  of a target region?"

### Practical considerations

- Document the climate data period used for TPE definition
- Test if TPE boundaries have shifted over recent decades
- Consider using climate projections (CMIP6) to forecast TPE shifts
- Report uncertainty: a TPE defined from 10 years of data is less robust
  than one from 30 years

## Species-Specific Notes

### Cereals (wheat, maize, rice, barley, sorghum)

- Most enviromics literature; well-characterized critical windows
- GDD base temperatures well established
- Photoperiod sensitivity varies widely among varieties of the same species
- Vernalization requirement adds complexity (winter wheat/barley)

### Legumes (soybean, common bean, chickpea, lentil)

- Photoperiod sensitivity often more critical than temperature
- Maturity group classification captures adaptation patterns
- Water availability during pod-fill is the dominant environmental
  constraint for yield

### Perennial crops (coffee, cacao, oil palm, fruit trees)

- Multi-year carryover effects: this year's environment affects next
  year's yield
- Age effects confound environmental effects
- Longer breeding cycles mean fewer MET data points

### Polyploid crops (wheat, potato, sugarcane, cotton)

- Subgenome-specific responses to environment are possible
- Homeologous gene expression may respond differently to stress
- Genomic selection with G×E is active research area

## Software Ecosystem

| Tool/Package | Purpose | License |
|-------------|---------|---------|
| ASReml-R | FA models, sparse MET | Commercial |
| Sommer | FA/Mixed models (open source) | GPL |
| metan | AMMI, GGE, stability, MTSI | GPL |
| agricolae | Classical stability, AMMI | GPL |
| lme4 | Mixed models, reaction norms | GPL |
| BreedR | Spatial + genetic models | GPL |
| EnvRtype | Environmental covariates, reaction norms | GPL |
| nasapower | NASA POWER weather data | MIT |
| geodata | WorldClim, SoilGrids access | GPL |
