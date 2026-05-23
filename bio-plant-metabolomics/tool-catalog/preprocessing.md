# Metabolomics Data Preprocessing

**Goal:** Convert raw vendor files to open formats (mzML), perform quality assessment, and prepare data for peak detection
**Best for:** All metabolomics workflows — essential first step before XCMS or MZmine

## Prerequisites
- ProteoWizard (msConvert) or ThermoRawFileParser
- R 4.0+ with xcms, MSnbase
- Raw MS files (.raw, .d, .wiff, .lcd)

## Raw File Conversion

### msConvert (ProteoWizard) — Universal Converter

```bash
# Convert Thermo .raw to mzML (centroided)
msconvert data.raw \
  --mzML \
  --filter "peakPicking true 1-" \
  --filter "msLevel 1" \
  --outdir mzML/

# Convert with MS2 (for annotation)
msconvert data.raw \
  --mzML \
  --filter "peakPicking true 1-" \
  --filter "msLevel 1-2" \
  --outdir mzML_ms2/

# Batch conversion
for f in *.raw; do
  msconvert "$f" --mzML --filter "peakPicking true 1-" --outdir mzML/
done
```

### Bruker .d Files

```bash
msconvert data.d --mzML --filter "peakPicking true 1-" --outdir mzML/
```

### AB Sciex .wiff Files

```bash
msconvert data.wiff --mzML --filter "peakPicking true 1-" --outdir mzML/
```

## Quality Assessment

### Check TIC (Total Ion Chromatogram)

```r
library(xcms)
library(MSnbase)

files <- list.files("mzML/", pattern = "*.mzML", full.names = TRUE)
raw_data <- readMSData(files, mode = "onDisk")

# Extract TIC
tic <- chromatogram(raw_data, aggregationFun = "sum")

# Plot TIC overlay
plot(tic, col = rep(c("red", "blue", "green"), length.out = length(files)))
```

### QC Sample Assessment

```r
# Calculate TIC variation in QC samples
qc_indices <- grep("QC", basename(files))
qc_tic <- tic[qc_indices]

# Coefficient of variation per feature should be < 30% in QCs
# (This is calculated after peak detection)
```

## Data Organization

Recommended directory structure:

```
project/
├── raw/                    # Original vendor files
├── mzML/                   # Converted mzML files
├── metadata.csv            # Sample metadata
├── scripts/
│   ├── 01_convert.R
│   ├── 02_peak_detection.R
│   └── 03_diff_analysis.R
└── results/
    ├── figures/
    └── tables/
```

### Metadata Format

```csv
sample_name,group,donor,tissue,time_point,run_order
QC1,QC,NA,leaf,NA,1
S01,Control,Plant1,leaf,T0,2
S02,Control,Plant2,leaf,T0,3
S03,Treatment,Plant3,leaf,T1,4
QC2,QC,NA,leaf,NA,5
...
```

## Peak Picking Parameters (Pre-Optimization)

### LC-MS (Q-TOF, Orbitrap)

| Parameter | HPLC | UPLC |
|-----------|------|------|
| ppm | 15-25 | 5-10 |
| peakwidth (s) | c(10, 60) | c(2, 20) |
| snthresh | 6-10 | 3-6 |
| noise | 500-1000 | 100-500 |

### GC-MS (Quadrupole, TOF)

| Parameter | GC-Q | GC-TOF |
|-----------|------|--------|
| m/z tolerance | 0.01 Da | 0.005 Da |
| peakwidth (s) | c(3, 30) | c(2, 15) |
| snthresh | 10-20 | 5-10 |

## Missing Value Patterns

After peak detection and alignment, categorize missing values:

```r
# Identify missing value patterns
missing_pattern <- is.na(feature_table)
missing_by_sample <- colSums(missing_pattern)
missing_by_feature <- rowSums(missing_pattern)

# Plot missing value distribution
barplot(sort(missing_by_feature), main = "Missing values per feature")
```

- **Missing in all replicates of one group** → likely real biological absence
- **Missing randomly across samples** → likely below detection limit
- **Missing in specific runs** → potential batch effect or injection failure

## Plant-Specific Considerations

- Plant tissue extracts often need dilution (1:10 to 1:100) before injection to avoid overloading
- For leaf samples, include a wash run between every 5-10 samples to clear column
- Pigment-rich extracts (spinach, kale): use C18 SPE cleanup before LC-MS
- Volatile compounds (terpenes in mint, basil): consider HS-SPME-GC-MS instead of solvent extraction
- Root exudates: collected in water/hydroponics; need lyophilization + reconstitution before analysis

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| msConvert not found | ProteoWizard not installed | Install via `brew install proteowizard` (macOS) or download from website |
| "Unsupported format" | Vendor lock-in or old format | Use vendor software to export mzML |
| Huge mzML files (>1GB each) | Profile mode data | Use peak picking during conversion |
| No peaks in some samples | Sample preparation or injection failure | Check TIC traces; re-run failed samples |
| Different RT between samples | Column aging or gradient drift | Use QC bracketing; run system suitability test |
