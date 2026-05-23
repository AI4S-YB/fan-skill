# Metabolite Annotation (SIRIUS + CSI:FingerID)

**Goal:** Putatively identify metabolite structures from MS/MS data using computational methods
**Best for:** LC-MS/MS data where reference standard comparison is not possible

## Prerequisites
- SIRIUS 5.0+ (CLI) and CSI:FingerID
- MS/MS data (mzML or mgf format)
- Feature table with m/z, RT, and MS/MS spectrum accessions

## SIRIUS Workflow

### 1. Export MS/MS Data

From XCMS/CAMERA:
```r
library(CAMERA)
# Annotate isotopes and adducts
xsa <- annotate(xdata, CAMERA_param)

# Export MS/MS spectra to MGF
export_ms2_spectra(xdata, xsa, file = "ms2_spectra.mgf")
```

From MZmine:
```
MZmine → Export → MS/MS spectra → MGF format
```

### 2. Run SIRIUS Molecular Formula Identification

```bash
# Process each MS/MS spectrum
sirius \
  --input ms2_spectra.mgf \
  --output sirius_results \
  --maxmz 800 \
  --elements "CHNOPS" \
  --ppm-max 10 \
  --candidates 10 \
  --processors 16 \
  --profile orbitrap   # or qtof
```

### 3. Run CSI:FingerID Structure Database Search

```bash
# Search against structure databases
sirius --input ms2_spectra.mgf \
  --output sirius_fingerid_results \
  --processors 16 \
  --fingerid \
  --databases "pubchem,bio,pubmed,plants" \
  --compound-cache sirius_results/compound_cache.sqlite
```

### 4. Run CANOPUS Compound Class Prediction

```bash
# Predict compound class (no database needed)
sirius --input ms2_spectra.mgf \
  --output sirius_canopus_results \
  --processors 16 \
  --canopus
```

## Interpreting SIRIUS Results

### Key Output Files

| File | Content |
|------|---------|
| formula_candidates.tsv | Ranked molecular formula candidates |
| structure_candidates.tsv | CSI:FingerID structure hits |
| canopus_summary.tsv | CANOPUS class predictions |
| compound_summary.tsv | Combined summary per feature |

### Confidence Metrics

```r
library(readr)

# Read SIRIUS results
formulas <- read_tsv("formula_candidates.tsv")
structures <- read_tsv("structure_candidates.tsv")
canopus <- read_tsv("canopus_summary.tsv")

# Filter for high-confidence annotations
confident <- structures %>%
  filter(
    CSI_FingerID_score > -100,    # Higher is better (closer to 0)
    confidence_score > 0.5,       # Overall confidence
    tanimoto_score > 0.3          # Structural similarity to DB hit
  )

# Best hit per feature
best_hits <- confident %>%
  group_by(feature_id) %>%
  slice_max(order_by = CSI_FingerID_score, n = 1)
```

## Alternative: GNPS Molecular Networking

```bash
# For community-based annotation via GNPS
# 1. Upload MGF to https://gnps.ucsd.edu/
# 2. Run METABOLOMICS-SNETS-V2 workflow
# 3. Retrieve network and library hits
```

## Manual Annotation Verification

For top candidate metabolites, verify against:
1. **Exact mass**: match within instrument ppm tolerance
2. **MS/MS fragments**: do they match reported fragments in literature?
3. **Retention time**: is the RT consistent with compound polarity?
4. **Biological context**: is this compound known in your plant species/family?
5. **Adduct pattern**: does the adduct(s) make chemical sense?

## Key Parameters

| Parameter | Purpose |
|-----------|---------|
| --ppm-max | Mass accuracy in ppm |
| --elements | Elements to consider in formula |
| --candidates | Number of formula candidates |
| --profile | Instrument type (orbitrap, qtof, ft-icr) |
| --databases | Structure databases for CSI:FingerID |

## Plant-Specific Databases

Supplement SIRIUS with plant-specific compound databases:
- **LOTUS** (https://lotus.naturalproducts.net/): Natural products database
- **COCONUT** (https://coconut.naturalproducts.net/): Collective open natural products database
- **CMAUP** (http://bidd.group/CMAUP/): Collective molecular activities of useful plants
- **NPASS** (http://bidd.group/NPASS/): Natural product activity & species source

## MSI Confidence Level Assignment

| Level | Evidence | How to Achieve |
|-------|----------|---------------|
| Level 1 | Standard match (RT + MS + MS/MS) | Compare with authentic standard |
| Level 2 | MS/MS spectral library match | CSI:FingerID score > -50 + manual review |
| Level 3 | Molecular formula + candidate | SIRIUS formula prediction + DB match |
| Level 4 | Molecular formula only | SIRIUS formula prediction, no structural hits |
| Level 5 | m/z + RT only | Feature detected but not annotated |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No MS/MS data" | DDA not acquired | Cannot annotate; only report Level 5 |
| SIRIUS timeout | Too many features | Filter to significant features only |
| Low CSI:FingerID scores | Compound not in database | Plant compound may be novel; use CANOPUS class prediction |
| Too many candidates per feature | Poor MS/MS quality | Increase collision energy; exclude low TIC spectra |
