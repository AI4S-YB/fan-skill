# fan-skill

AI-powered plant bioinformatics and breeding analysis engine. From biological question to publication-ready results. Three-layer knowledge architecture (personal → general → Agent fallback) with B+C dual-mode decision system.

> 🇨🇳 [中文文档](README_CN.md)

## Quick Start

```bash
git clone git@github.com:AI4S-YB/fan-skill.git
cd fan-skill

# Install for your AI coding agent:
bash install-claude.sh     # Claude Code
bash install-codex.sh      # Codex CLI
bash install-gemini.sh     # Gemini CLI
bash install-opencode.sh   # OpenCode
bash install-all.sh        # Auto-detect and install for all

# Or use the universal installer:
# npx skills add AI4S-YB/fan-skill
```

Then describe your biological question:

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
│  Intent → Three-Layer Match → Execute
└──────────────┬──────────────────┘
               │
    ┌──────────┼──────────┐
    ▼          ▼          ▼
┌────────┐ ┌───────┐ ┌───────┐
│  User  │ │ gwas  │ │rnaseq │  ← user-knowledge/ (your experience)
│ entry  │ │ rules │ │ rules │     knowledge-base/ (31 curated)
└───┬────┘ └───┬───┘ └───┬───┘
    │          │         │
    └──────────┼─────────┘
               ▼
┌─────────────────────────────────┐
│  engine/                        │  ← Shared execution engine
│  rule_engine · validate · log   │
└─────────────────────────────────┘
               │
    ┌──────────┼──────────┐
    ▼          ▼          ▼
  分析报告   决策日志   图表
               │
    (未匹配时) Agent 兜底推理 → 沉淀到 user-knowledge/
```

Validated on apple RNA-seq and rice candidate gene analysis. Each decision is auditable.

### Three-Layer Knowledge

| Layer | Path | Priority | Content |
|-------|------|:--------:|---------|
| **User** | `user-knowledge/` | 🔺 Highest | Your personal analysis experience. Starts empty, grows with use. |
| **General** | `knowledge-base/` | 🔸 Default | 31 curated entries covering major plant bioinfo workflows. |
| **Agent Fallback** | (model reasoning) | 🔻 Last resort | When neither layer matches, the Agent reasons and can precipitate results into your user layer. |

### B+C Dual-Mode

| Layer | Role | Form |
|------|------|------|
| **C-layer (Rule)** | Deterministic — same data → same method | `rules.yaml` |
| **B-layer (Expert)** | Flexible — expert reasoning for edge cases | `notebook.md` |
| **User control** | `decision_mode: rule | expert | hybrid` | `params.yaml` |
| **User entries** | Minimum `meta.yaml` + `notebook.md`; optional `rules.yaml` | `user-knowledge/` |

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
- **Fan-skill solution**: Three-layer knowledge architecture. Your personal experience (`user-knowledge/`) takes priority over the 31 curated entries (`knowledge-base/`). When neither matches, Agent reasoning serves as fallback and can be saved as new personal entries — the system grows with you.
- **B+C dual-mode**: C-layer encodes "what method" as structured rules; B-layer encodes "why" as expert reasoning. Both layers coexist in every entry.

Every decision is auditable — `engine/log_decision.sh` records why each choice was made.

## Contribute

**To the general knowledge base:** New analysis = 4 files (not 14). See **[Contributing Guide](docs/en/contributing.md)**.

```bash
cp templates/rules-template.yaml knowledge-base/<name>/rules.yaml
cp templates/notebook-template.md knowledge-base/<name>/notebook.md
# + consult-guide.md + analysis-primer.md
engine/validate_entry.sh knowledge-base/<name>/
```

**To your personal knowledge base:** Just use fan-skill. When Agent fallback solves a novel analysis, choose to save it — an entry is auto-generated in `user-knowledge/drafts/`. Review and confirm to move it to `user-knowledge/confirmed/`.

## License

MIT
