# Peak Detection and Alignment (XCMS)

**Goal:** Detect chromatographic peaks and align features across LC-MS or GC-MS runs
**Best for:** Untargeted metabolomics — XCMS for LC-MS, MZmine for GC-MS

## Prerequisites
- R 4.0+ with xcms, MSnbase, CAMERA packages
- Converted mzML files (see `tool-catalog/preprocessing.md`)
- Sample metadata (phenotype/condition labels)

## XCMS Workflow (LC-MS)

### 1. Load Data and Define Phenotype

```r
library(xcms)
library(MSnbase)

# Read mzML files
files <- list.files("mzML/", pattern = "*.mzML", full.names = TRUE)
pd <- data.frame(
  sample_name = sub(".mzML", "", basename(files)),
  group = c(rep("Control", 5), rep("Treatment", 5), rep("QC", 3))
)
raw_data <- readMSData(files, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")
```

### 2. Peak Detection (CentWave)

```r
# Define parameters
cwp <- CentWaveParam(
  ppm = 15,               # Mass tolerance (adjust for instrument: Q-TOF=25, Orbitrap=5)
  peakwidth = c(5, 20),   # Peak width range in seconds (UPLC: 2-10, HPLC: 5-30)
  snthresh = 5,           # Signal-to-noise threshold
  prefilter = c(3, 100),  # At least 3 consecutive scans > 100 intensity
  mzCenterFun = "wMean",
  noise = 500
)

# Run peak detection
xdata <- findChromPeaks(raw_data, param = cwp)
```

### 3. Retention Time Alignment (Obiwarp)

```r
# Group peaks across samples
pgp <- PeakDensityParam(
  sampleGroups = xdata$group,
  bw = 5,                 # Bandwidth for grouping
  minFraction = 0.5       # Feature must be present in >50% of samples in at least one group
)
xdata <- groupChromPeaks(xdata, param = pgp)

# Retention time correction
obp <- ObiwarpParam(
  binSize = 0.6,
  centerSample = which(xdata$group == "QC")[1]  # Use QC as reference
)
xdata <- adjustRtime(xdata, param = obp)

# Re-group after alignment
xdata <- groupChromPeaks(xdata, param = pgp)
```

### 4. Fill Missing Peaks

```r
# Fill peaks that were missed in some samples
xdata <- fillChromPeaks(xdata)
```

### 5. Extract Feature Table

```r
feature_table <- featureValues(xdata, value = "into")  # Peak area
feature_info <- featureDefinitions(xdata)               # m/z and RT

# Log2 transform
feature_table_log2 <- log2(feature_table + 1)
write.csv(feature_table_log2, "feature_table_log2.csv")
write.csv(feature_info, "feature_info.csv")
```

## MZmine Workflow (GC-MS)

For GC-MS data, MZmine is preferred due to its built-in EI spectral library search:

### Graphical Workflow (MZmine 3)

1. **Import**: Raw data → mzML files
2. **Mass Detection**: Centroid → Noise level 1.0E3 (GC-MS)
3. **Chromatogram Builder**: Min time span 0.1min, Min height 1.0E3
4. **Deconvolution**: Local minimum search or ADAP
5. **Isotope Grouping**: Remove isotopes
6. **Alignment**: Join Aligner (m/z tolerance 0.01 Da, RT tolerance 0.1 min)
7. **Gap Filling**: Peak finder (intensity tolerance 20%)
8. **Export**: CSV feature table

### Batch Mode (MZmine CLI)

```bash
mzmine -batch batch_workflow.xml
```

## Key Parameters Comparison: XCMS vs MZmine

| Parameter | XCMS (CentWave) | MZmine (ADAP) |
|-----------|-----------------|----------------|
| Mass tolerance | ppm (instrument-dependent) | m/z tolerance (Da) |
| Peak width | c(min, max) seconds | Min/max peak duration (min) |
| S/N threshold | snthresh | Noise level |
| Alignment method | Obiwarp (profile-based) | Join Aligner (match-based) |

## Plant-Specific Considerations

- Plant secondary metabolites produce many isomeric peaks — use tight RT tolerance
- Chlorophyll derivatives produce intense background in leaf extracts — may need to exclude known m/z regions
- Plant phenolics often form adducts ([M-H]-, [M+HCOO]-, [2M-H]-) — use CAMERA for adduct annotation
- Glycosylated metabolites (flavonoid glycosides) show characteristic neutral losses in MS/MS

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "No peaks detected" | snthresh or peakwidth too strict | Relax snthresh; check raw data quality |
| "Error in align: too few peaks" | Insufficient common features | Lower minFraction; check if samples are from same platform |
| Excessive features (>10,000) | Noise picked up as peaks | Increase snthresh or prefilter steps |
| RT alignment failed | Extreme chromatographic drift | Check LC column condition; exclude outlier runs |
| Peak area all zeros after fill | FillChromPeaks integration failed | Increase integration window; check raw data |
