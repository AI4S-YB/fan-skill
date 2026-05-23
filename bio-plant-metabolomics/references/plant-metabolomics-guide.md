# Plant Metabolomics Guide

## Why Plant Metabolomics Is Distinct

Plant metabolomics differs fundamentally from clinical/metabolic metabolomics:

| Aspect | Clinical Metabolomics | Plant Metabolomics |
|--------|----------------------|-------------------|
| Metabolome size | ~3,000-5,000 endogenous | 20,000-200,000+ (species-dependent) |
| Known metabolites | > 80% in databases | < 20% in databases |
| Chemical diversity | Moderate | Extreme (alkaloids, flavonoids, terpenoids, etc.) |
| Reference standards | Widely available | Limited (many are species-specific) |
| Biological replicates | Tightly controlled (lab animals/cell lines) | High variability (field-grown plants) |
| Temporal dynamics | Relatively stable | Strong diurnal and developmental changes |
| Matrix effects | Moderate | High (pigments, polyphenols, lipids) |

## Platform Selection Guide

### Choose LC-MS When:
- Interested in secondary/specialized metabolites (flavonoids, alkaloids, glucosinolates, terpenoids)
- Broad metabolome coverage needed
- High sensitivity required
- Have access to MS/MS for structural characterization

### Choose GC-MS When:
- Focused on primary metabolism (sugars, organic acids, amino acids, fatty acids)
- Want robust quantification
- Plant volatile analysis (headspace GC-MS)
- Need high reproducibility and established spectral libraries

### Choose NMR When:
- Absolute quantification needed
- Non-destructive analysis (e.g., in vivo metabolomics)
- Structural elucidation of novel compounds
- Small number of expected metabolites (fingerprinting)

### Complementary Approaches (Recommended)
Run both LC-MS and GC-MS on the same extracts for comprehensive coverage:
- GC-MS: captures primary metabolism
- LC-MS: captures secondary/specialized metabolism
- Combined: covers > 90% of detectable metabolome

## Extraction Methods

### General Polar Extraction (Methanol:Water, 80:20)
- Good for flavonoids, phenolic acids, amino acids, sugars
- Standard for untargeted plant metabolomics
- 50-100mg fresh weight tissue in 1mL solvent
- Homogenize (bead beater) + sonicate (15 min) + centrifuge

### Non-polar Extraction (Chloroform:Methanol:Water)
- Good for lipids, chlorophyll, carotenoids
- Phase separation yields polar + non-polar fractions

### Volatile Extraction
- HS-SPME (headspace solid-phase microextraction) for terpenes
- Direct thermal desorption for volatile organic compounds (VOCs)

### Acidified Extraction
- 0.1% formic acid in extraction solvent
- Stabilizes anthocyanins (flavylium cation form)
- Improves extraction of alkaloids (protonation)

## Replication Requirements

Plant metabolomics needs MORE replicates than mammalian metabolomics:

| Study Type | Minimum Replicates | Recommended |
|------------|-------------------|-------------|
| Field-grown plants | 5-6 | 8-10 |
| Greenhouse plants | 4-5 | 6-8 |
| Growth chamber | 3-4 | 5-6 |
| Cell/tissue culture | 3 | 4-5 |

Rationale: Environmental variation (light, humidity, soil heterogeneity) adds substantial metabolome variance in plants.

## Quality Control (QC) Strategy

### Pooled QC Samples
- Mix equal aliquots of ALL study samples
- Inject at beginning (5-10 conditioning runs), every 5-10 study samples, and at end
- CV of features in QCs should be < 30%

### Internal Standards
- Add known amount of non-endogenous compound to each sample before extraction
- Examples: lidocaine (ESI+), camphor sulfonic acid (ESI-), ribitol (GC-MS)
- Monitor IS signal across batch for drift

### Blank Extraction
- Process solvent-only "sample" through entire extraction pipeline
- Identifies background signals from solvents, tubes, and columns
- Features common in blanks should be filtered

## Metabolite Identification Strategy for Plants

### Tier 1: Database-First Approach
1. Match accurate mass against plant-specific databases (LOTUS, COCONUT, PlantCyc)
2. Filter by taxonomy (is this compound reported in your plant family?)
3. Check MS/MS fragmentation against spectral libraries (MoNA, GNPS, MassBank)
4. Verify with literature (has this compound been reported in your species?)

### Tier 2: Computational Annotation
1. SIRIUS: predict molecular formula from isotope pattern + MS/MS
2. CSI:FingerID: search structure databases with predicted fingerprints
3. CANOPUS: predict compound class (even without database match)

### Tier 3: Manual Interpretation
1. Analyze MS/MS fragmentation pattern for characteristic neutral losses
2. Compare with known analogs
3. Consider biosynthetic logic (e.g., flavonoids typically glycosylated at 3-O or 7-O positions)

## Species-Specific Metabolic Features

### Brassicaceae (Arabidopsis, Brassica, Raphanus)
- Glucosinolates and hydrolysis products (isothiocyanates, nitriles)
- Hydroxycinnamic acid derivatives (sinapic acid esters)

### Fabaceae (Soybean, Medicago, Phaseolus)
- Isoflavonoids (daidzein, genistein, coumestrol)
- Saponins (triterpenoid and steroidal)

### Solanaceae (Tomato, Potato, Pepper, Tobacco)
- Steroidal glycoalkaloids (tomatine, solanine, chaconine)
- Acyl sugars (glandular trichome metabolites)

### Lamiaceae (Mint, Basil, Rosemary, Lavender)
- Volatile terpenoids (menthol, thymol, carvacrol)
- Rosmarinic acid and derivatives

### Poaceae (Rice, Maize, Wheat, Barley)
- Benzoxazinoids (DIMBOA, MBOA) — defense compounds
- Phenolic amides (avenanthramides in oat)
