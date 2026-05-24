# Fan-Skill Quick Start / 快速开始

## Installation / 安装

```bash
git clone git@github.com:AI4S-YB/fan-skill.git
cd fan-skill
bash install.sh
```

## First Use / 第一次使用

Open Claude Code and describe your biological question in natural language.
在 Claude Code 中直接用自然语言描述你的生物学问题：

> "I have 300 rice accessions with GBS data and 3 years of yield. Find genes controlling grain weight."
> "我有 300 份水稻材料的 GBS 数据和 3 年产量表型，想找到控制粒重的基因。"

Fan-skill will automatically match your intent to the knowledge base, check your
experimental design, recommend an analysis path, and execute it.
Fan-skill 会自动匹配知识库、检查实验设计、推荐分析路径并执行。

## Usage Modes / 使用方式

| Mode / 模式 | When / 场景 | Example / 示例 |
|------|------|------|
| **I know what analysis I need** / 我知道要做什么 | You have data and a clear analysis goal | "Run GWAS on my VCF and phenotype data" |
| **I have data, what can I do?** / 我有数据，能做什么？ | You have data but aren't sure what analyses are possible | "I have RNA-seq counts from drought-treated rice. What can I analyze?" |
| **I have a biological question** / 我有生物学问题 | You want to answer a research question but don't know what analysis to use | "Find genes controlling salt tolerance in my breeding population" |

## Requirements / 环境要求

- Claude Code (any recent version)
- R 4.0+ (for most analyses / 大部分分析需要)
- Python 3.10+ (for data inspection and GFF3 parsing / 数据探查和 GFF3 解析)
- Individual analyses may require additional tools (PLINK, GAPIT, BWA, etc.).
  Run `engine/check_env.sh` to see what's available.
  个别分析需要额外工具。运行 `engine/check_env.sh` 查看可用工具。

## Supported Species / 支持的物种

Rice 水稻, Maize 玉米, Wheat 小麦, Soybean 大豆, Cotton 棉花, Rapeseed 油菜,
Arabidopsis 拟南芥, Tomato 番茄, Potato 马铃薯, Barley 大麦, Sorghum 高粱,
Sugarcane 甘蔗 — and any species with a reference genome.
以及任何有参考基因组的物种。

## What's Next / 下一步

- See all 30 analysis capabilities: `docs/capability-catalog.md`
- Contribute a new analysis: `docs/contributing.md`
- Understand the architecture: `docs/superpowers/specs/2026-05-22-fan-skill-design.md`
