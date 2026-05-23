# QTL IciMapping — Inclusive Composite Interval Mapping (ICIM) & MET QTL

**Goal:** Integrated QTL mapping software — genetic map construction + ICIM + multi-environment QTL analysis
**Best for:** High-density biparental populations, multi-environment trials (MET), QTLxE analysis

## Prerequisites
- QTL IciMapping 4.2+ (standalone software, Windows/Linux)
- Genotype: tab-delimited, markers x individuals
- Phenotype: tab-delimited, individuals x traits, multiple environments supported
- Map file: chromosome and position for each marker (or let software build map)

## Basic Usage

### Input File Format

Genotype file (`.bip` format for biparental):
```
# Header: marker names
# Rows: individuals (first column = individual ID)
# Values: 0 (homozygous parent1), 1 (heterozygous), 2 (homozygous parent2), -1 (missing)
```

Phenotype file (`.txt`):
```
# Header: trait names (if multi-env: trait_env1, trait_env2, ...)
# Rows: individuals
```

### Workflow in IciMapping GUI

1. **Create project** → Select population type (F2/RIL/DH/BC)
2. **Import data** → Genotype file + phenotype file
3. **Build genetic map** (optional if map file provided):
   - `Grouping` → Assign markers to linkage groups (LOD threshold usually 3-10)
   - `Ordering` → Order markers within linkage groups (nnTwoOpt, SER, or RECORD)
   - `Rippling` → Refine marker order
4. **ICIM mapping**:
   - Select mapping method: ICIM-ADD (additive effect) or ICIM-DOM (both additive + dominance)
   - Set scanning step (1 cM recommended)
   - Set PIN (probability in stepwise regression): 0.001 for inclusion, 0.002 for exclusion
   - Run scanning
5. **MET QTL analysis** (if multi-environment):
   - Select `MET` option
   - Choose joint or per-environment output
   - Estimate QTL main effect + QTLxE interaction

### Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| LOD grouping threshold | 3-5 (low-density), 10+ (high-density) | Higher = fewer false linkage groups |
| Mapping function | Kosambi (standard) or Haldane (no interference) | Kosambi more realistic for plants |
| Step size | 1 cM | Standard for QTL scanning |
| PIN (ICIM) | 0.001 | Default inclusion threshold is appropriate |
| LOD threshold (MET) | 2.5-3.0 | Slightly lower than single-env due to increased df |

## ICIM vs CIM (R/qtl)

| Aspect | ICIM (IciMapping) | CIM (R/qtl) |
|--------|-------------------|-------------|
| Background control | Marker selection via stepwise regression | User-specified cofactors |
| Multi-env support | Native MET module | Requires separate scripts |
| Ease of use | GUI-based, beginner-friendly | R scripting, more flexible |
| Statistical power | Slightly higher for linked QTLs | Comparable for isolated QTLs |
| Publication | Common in crop journals | More common in genetics journals |

## MET QTL Analysis

Two strategies in IciMapping:

### Strategy 1: Joint Mapping
- Estimates QTL main effect (across environments)
- Estimate QTLxE variance for each QTL
- More statistical power than per-environment scanning

### Strategy 2: Per-Environment Mapping
- Run ICIM independently in each environment
- Compare QTL positions across environments
- Identify stable vs environment-specific QTL

## Plant-Specific Notes

- **ICIM is the default in Chinese crop QTL studies** — if publishing in a Chinese journal, IciMapping is often the expected tool.
- **MET module**: Excellent for rice, wheat, maize, soybean MET data where the same population was phenotyped in multiple locations.
- **Map construction within IciMapping**: Reasonable for small-medium datasets but slower than R/qtl est.map. For >5000 markers, build map externally and import.
- **Output interpretation**: ICIM reports LOD peak position + left/right interval boundaries, which map to physical positions if a reference genome is available.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Cannot read file" | Wrong file format | Ensure tab-delimited, no empty header cells |
| "No linkage groups found" | LOD threshold too high | Lower grouping LOD to 3-5 |
| "Negative distances" | Markers in wrong order | Use rippling function or check input map |
| "Overly long chromosomes" | Inflated recombination distances | Check marker order; remove distorted markers |
| "MET analysis fails" | Missing environments in phenotype | Ensure all individuals have phenotype for all envs |
