# Fan-Skill Quick Start

## Installation

```bash
git clone git@github.com:AI4S-YB/fan-skill.git
cd fan-skill
bash install.sh
```

## First Use

Open Claude Code and describe your biological question in natural language:

> "I have 300 rice accessions with GBS data and 3 years of yield. Find genes controlling grain weight."

Fan-skill will automatically match your intent to the knowledge base, check your experimental design, recommend an analysis path, and execute it.

## Usage Modes

| Mode | When | Example |
|------|------|------|
| **I know what I need** | You have data and a clear analysis goal | "Run GWAS on my VCF and phenotype data" |
| **I have data, what can I do?** | You have data but aren't sure what analyses are possible | "I have RNA-seq counts from drought-treated rice" |
| **I have a biological question** | You want to answer a research question but don't know the analysis | "Find genes controlling salt tolerance" |

## Requirements

- Claude Code (any recent version)
- R 4.0+ (for most analyses)
- Python 3.10+ (for data inspection and GFF3 parsing)
- Individual analyses may require additional tools (PLINK, GAPIT, BWA, etc.). Run `engine/check_env.sh` to check.

## Supported Species

Rice, maize, wheat, soybean, cotton, rapeseed, Arabidopsis, tomato, potato, barley, sorghum, sugarcane — and any species with a reference genome.

## What's Next

- See all 30 analysis capabilities: [Capability Catalog](capability-catalog.md)
- Contribute a new analysis: [Contributing Guide](contributing.md)
- Understand the architecture: [Design Spec](../superpowers/specs/2026-05-22-fan-skill-design.md)
