# Contributing to Fan-Skill

Thank you for contributing to fan-skill. This guide explains how the project is structured and how to add a new analysis capability.

## Architecture: B+C Pattern

Fan-skill uses a **B+C (Base + Configuration)** architecture. Every analysis entry consists of two layers:

- **B (Base)** — a ready-to-run analysis notebook (`.md`) containing complete R/Python code, plant-specific design gates, QC checkpoints, and troubleshooting hooks. This is what gets executed when you invoke the skill.
- **C (Configuration)** — a `rules.yaml` file that describes when and how to use the analysis, including intent matching, input validation, design gate logic, and integration with other entries.

This separation means that subject-matter experts write the analysis code (B), while the skill's inference engine uses the configuration (C) to select, customize, and execute the right analysis for the user's context.

## Project Structure

```
fan-skill/
├── knowledge-base/          # B+C entries, one directory per analysis
│   ├── gwas/                # Example entry
│   │   ├── rules.yaml       # C: configuration
│   │   └── notebook.md      # B: analysis notebook
│   ├── rnaseq/
│   ├── genomic-selection/
│   └── ... (27+ more)
├── tool-registry/           # Tool wrappers (PLINK, GAPIT, BWA, etc.)
├── references/             # Reference data (species cheatsheet, QC thresholds)
├── templates/              # Reusable document templates
├── engine/                 # Core engine scripts (intent matching, env checks)
├── docs/                   # Documentation (you are here)
│   ├── en/                 # English docs
│   └── zh/                 # Chinese docs
└── theme/                  # Output styling
```

## Adding a New Analysis: 4 Files

To contribute a new analysis capability, create these 4 files:

### 1. `rules.yaml` — Configuration

```yaml
# knowledge-base/<entry-name>/rules.yaml
entry: <entry-name>
display_name: <Human-readable name>
category: <one of the capability catalog categories>

intent:
  triggers:             # Keywords/phrases that trigger this entry
    - <trigger phrase 1>
    - <trigger phrase 2>
  anti_triggers:        # Keywords that suppress this entry
    - <anti-trigger>

inputs:
  required:
    - name: <input-name>
      type: <file-type>
      description: <what this is>
      format: [<accepted formats>]
  optional:
    - name: <input-name>
      type: <file-type>
      description: <what this is>

design_gates:           # Plant-specific experimental design checks
  - gate: <gate-name>
    check: <what to verify>
    action: <what to do if check fails>

outputs:
  - <output description 1>
  - <output description 2>

depends_on:             # Optional: other entries this one depends on
  - <entry-name>

related:                # Optional: related entries
  - <entry-name>
```

### 2. `notebook.md` — Analysis Notebook

The notebook is a literate programming document that combines explanation with executable code. It must include:

- **Objective**: What this analysis does and when to use it
- **Plant-specific design gates**: Check ploidy, breeding system, population structure before running
- **Input checklist**: Validate all inputs with format checks
- **Step-by-step analysis**: R or Python code blocks with inline explanation
- **QC checkpoints**: At key steps, verify intermediate results (e.g., MAF distribution, PCA outliers)
- **Output tables and figures**: Publication-quality visualizations and formatted tables
- **Troubleshooting hooks**: Common failure modes and how to recover

Code blocks use the standard fenced format:

````markdown
```r
# R code here
library(ggplot2)
...
```
````

### 3. `tool-registry/<tool-name>.md` — Tool Wrapper (if needed)

If your analysis uses a tool not already in the tool registry, create a wrapper:

```markdown
# <Tool Name>

## Installation
...

## Parameter Reference
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|

## Parameter Decision Table
| Scenario | Recommended Settings | Rationale |
|----------|---------------------|-----------|

## Usage
...
```

### 4. Update `capability-catalog.md`

Add your entry to the appropriate category table in both `docs/en/capability-catalog.md` and `docs/zh/capability-catalog.md`.

## Design Gate Pattern

Design gates are a critical feature of fan-skill. They encode plant-specific experimental design knowledge and prevent inappropriate analyses. Every notebook must implement design gates at the top.

Example: before running GWAS on a crop species, check:

1. **Ploidy gate**: Is this an auto-polyploid? If yes, recode genotypes as diploid or use polyploid-aware methods.
2. **Population structure gate**: Is there uncontrolled structure? If yes, include PCA covariates.
3. **Self/cross-pollinated gate**: Self-pollinated → use FarmCPU/BLINK. Cross-pollinated → consider MLM with kinship.
4. **LD decay gate**: Estimate LD decay distance for significance threshold adjustment.

Each gate has three components in `rules.yaml`:

```yaml
design_gates:
  - gate: <gate-name>
    check: <what to verify>
    action: <what to do if check fails>
```

And a corresponding implementation in `notebook.md`:

```r
# ---- Design Gate: Population Structure ----
# Check: Is there uncontrolled population structure?
pca <- read.table("pca.eigenvec", header=FALSE)
var_explained <- pca_eigenval / sum(pca_eigenval) * 100
if (var_explained[1] > 10) {
  message("WARNING: PC1 explains >10% variance.")
  message("ACTION: Including top PCs as covariates.")
  covariates <- pca[, 1:3]  # Use top 3 PCs
}
```

## Pull Request Process

1. **Fork** the repository and create a feature branch
2. **Create the 4 files** described above for your new analysis
3. **Validate** your `rules.yaml` against the schema:
   ```bash
   python engine/validate_rules.py knowledge-base/<entry-name>/rules.yaml
   ```
4. **Test** your notebook by running it end-to-end with the example data provided in `test_matrix.csv`
5. **Check for integration issues**: your entry should not conflict with existing entries' intent triggers
6. **Submit a PR** with a description of:
   - What analysis you are adding
   - What plant-specific design gates you implemented
   - Any new tool-registry entries
   - Test results

## Style Guidelines

- **Notebook code**: Use `set.seed()` for reproducibility. Comment in English for R, English or Chinese for markdown prose.
- **rules.yaml**: All keys and values in English (the inference engine parses this programmatically).
- **Design gates**: Gates should be constructive (suggest a fix), not blocking (reject with no alternative).
- **References**: Cite methods papers where appropriate. Use standard journal abbreviations.

## Questions?

Open an issue on GitHub or start a discussion. We are especially interested in contributions for non-model plant species, polyploid-specific methods, and integration with public plant databases (Ensembl Plants, Phytozome, Gramene, etc.).
