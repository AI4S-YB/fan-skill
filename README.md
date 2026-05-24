# fan-skill

AI-powered plant bioinformatics and breeding analysis engine. From biological question to publication-ready results. Knowledge-base driven with B+C dual-mode decision architecture.

> 🇨🇳 [中文文档](README_CN.md)

## Quick Start

```bash
git clone git@github.com:AI4S-YB/fan-skill.git
cd fan-skill
bash install.sh
```

Open Claude Code and describe your biological question:

> "I have 300 rice accessions with GBS data and 3 years of yield. Find genes controlling grain weight."

Fan-skill will match your intent to the knowledge base, check your experimental design, recommend an analysis path, and execute it.

📖 **[Quick Start Guide](docs/en/quick-start.md)** · **[Analysis Capabilities](docs/en/capability-catalog.md)** · **[Contributing](docs/en/contributing.md)**

## Architecture

```
User says "find genes controlling grain weight"
        │
        ▼
┌─────────────────────────────────┐
│  SKILL.md                       │  ← Unified entry point
│  Intent → Match → Discover → Execute
└──────────────┬──────────────────┘
               │
    ┌──────────┼──────────┐
    ▼          ▼          ▼
┌───────┐ ┌───────┐ ┌───────┐
│ gwas  │ │rnaseq │ │  ...  │  ← knowledge-base/ (30 entries)
│ rules │ │ rules │ │       │     Each: C-layer rules + B-layer notebook
└───┬───┘ └───┬───┘ └───┬───┘
    │         │         │
    └─────────┼─────────┘
              ▼
┌─────────────────────────────────┐
│  engine/                        │  ← Shared execution engine
│  rule_engine · validate · log   │
└─────────────────────────────────┘
              │
              ▼
    Analysis report + decision log + figures
```

Validated on apple RNA-seq and rice candidate gene analysis. Each decision is auditable.

### B+C Dual-Mode

| Layer | Role | Form |
|------|------|------|
| **C-layer (Rule)** | Deterministic — same data → same method | `rules.yaml` |
| **B-layer (Expert)** | Flexible — expert reasoning for edge cases | `notebook.md` |
| **User control** | `decision_mode: rule | expert | hybrid` | `params.yaml` |

## Capabilities

| Goal | Fan-skill can |
|------|-------------|
| Find trait genes | GWAS · QTL mapping · eQTL · candidate gene association |
| Understand function | RNA-seq DE · time-series · GRN · multi-omics |
| Predict breeding value | Genomic selection · hybrid prediction · enviromics |
| Process raw data | Variant calling · genotype imputation · phenotype analysis |
| Epigenomics | ChIP-seq · ATAC-seq · DNA methylation |
| Microbiome | 16S/ITS amplicon · metagenomics |
| Metabolomics & Proteomics | LC-MS/GC-MS · DIA/DDA · PPI |
| Evolution | Comparative genomics · pan-genome |
| Edit & breed | CRISPR design · marker development |
| Visualize | Publication-quality ggplot2 figures |

Full catalog: **[Analysis Capabilities](docs/en/capability-catalog.md)**

## Design

Fan-skill encodes **judgment knowledge** — how an expert decides what analysis to run and which method to use. Not fixed pipelines.

- **Pipeline problem**: fixed steps fail on diverse data
- **Atomic skill problem**: free composition creates unstable invocation paths
- **Fan-skill solution**: C-layer encodes "what method" as structured rules; B-layer encodes "why" as expert reasoning. Users control stability vs flexibility via `decision_mode`.

Every decision is auditable — `engine/log_decision.sh` records why each choice was made.

## Contribute

New analysis = 4 files (not 14). See **[Contributing Guide](docs/en/contributing.md)**.

```bash
cp templates/rules-template.yaml knowledge-base/<name>/rules.yaml
cp templates/notebook-template.md knowledge-base/<name>/notebook.md
# + consult-guide.md + analysis-primer.md
engine/validate_entry.sh knowledge-base/<name>/
```

## License

MIT
