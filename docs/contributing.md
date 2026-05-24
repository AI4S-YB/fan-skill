# Fan-Skill Contributing Guide / 贡献指南

## Adding a New Analysis / 新增分析

Each analysis needs 4 files. Start from templates:
每个分析需要 4 个文件。从模板开始：

```bash
cp templates/rules-template.yaml knowledge-base/<your-analysis>/rules.yaml
cp templates/notebook-template.md knowledge-base/<your-analysis>/notebook.md
# Then create:
#   knowledge-base/<your-analysis>/consult-guide.md
#   knowledge-base/<your-analysis>/analysis-primer.md
```

### File Checklist / 文件清单

| File / 文件 | Minimum / 最低要求 | Purpose / 用途 |
|------|:--:|------|
| `consult-guide.md` | ≥30 lines / 行 | When to recommend this analysis, what to ask users / 何时推荐、问什么 |
| `rules.yaml` | ≥3 decision nodes / 决策节点 | C-layer: method selection rules + design gates / 方法选择规则+门禁 |
| `notebook.md` | ≥200 lines / 行 | B-layer: expert reasoning, pitfalls, plant-specific depth / 专家推理+陷阱 |
| `analysis-primer.md` | ≥20 lines / 行 | Plain-language result interpretation / 通俗结果解读 |

### rules.yaml Format / rules.yaml 格式

```yaml
name: <short-name>         # e.g. gwas, rnaseq
description: "<one-line>"  # What this analysis does
triggers: ["<keyword1>", "<keyword2>", ...]  # Match user intent
inputs: [<required_data>]  # What data is needed
outputs: [<produced_files>] # What this analysis produces

design_gates:              # Experimental design constraints
  - id: "<gate-id>"
    rule: "<requirement>"
    check: "<condition>"
    severity: "block | warn"
    on_fail: "<plain language explanation>"

<decision_node>:           # At least 3 decision nodes
  - rule_id: "<entry>-<method>-<number>"
    priority: 10           # Higher = evaluated first
    condition:
      <field>: "<value>"
    recommend: "<method>"
    reason: "<why>"
    tool_id: "<tool-name>" # Reference to tool-registry/

  - rule_id: "<entry>-fallback-999"
    priority: 0
    condition: "none_matched"
    action: "delegate_to_expert"
```

### notebook.md Quality Standards / notebook.md 质量标准

Notebooks guide the Agent's B-layer reasoning. A good notebook:
好的 notebook 引导 Agent 的 B 层推理。好的 notebook：

- Starts with the analysis mindset, not tool commands / 从分析思维开始，不是工具命令
- Uses plant-specific context (species, ploidy, breeding system) / 使用植物特化上下文
- Covers common pitfalls with actionable remedies / 覆盖常见陷阱并给出可操作的补救
- References related entries for chain analysis / 引用相关条目做链式分析
- Written in the language of the target users (Chinese for Chinese-speaking plant scientists) / 用目标用户的语言写作

### Tool Documentation / 工具文档

New tools used by your analysis should be documented in `tool-registry/`. Three levels:
新工具应记录在 `tool-registry/`。三级标准：

| Level / 级别 | Content / 内容 | When / 何时 |
|------|------|------|
| **Cookbook** | Code skeleton + parameter decisions + plant notes + errors | Core tools (DESeq2, GAPIT, PLINK...) |
| **Parameter Decision** | Existing params + "when to change" column | Secondary tools |
| **Basic** | Goal + basic usage + common errors | CLI wrappers, simple tools |

## Validation / 验证

Before submitting a PR, run / 提交 PR 前运行：

```bash
engine/validate_entry.sh knowledge-base/<your-analysis>/
```

This checks: 4 required files present, line counts, rules.yaml structure, rule_id on every rule.
检查: 4 文件齐全、行数达标、rules.yaml 结构完整、每条规则有 rule_id。

## Architecture Principles / 架构原则

Fan-skill encodes **judgment knowledge** — how an expert thinks about analysis decisions.
Fan-skill 编码的是**判断知识**——专家如何思考分析决策。

- **C-layer (rules.yaml)**: Structured, deterministic, testable. Same data → same method.
  C 层: 结构化、确定性的、可测试的。同样的数据 → 同样的方法。
- **B-layer (notebook.md)**: Narrative, flexible, expert reasoning. Handles edge cases.
  B 层: 叙事式、灵活的、专家推理。处理规则盲区。
- **B-C synergy**: When B and C diverge, B takes priority — but the override must be logged.
  B-C 协同: B 和 C 分歧时，B 优先——但覆盖必须记录。

The architecture is validated by 2 real-world tests (apple RNA-seq, rice candidate gene analysis).
架构经过 2 次真实测试验证（苹果 RNA-seq、水稻候选基因分析）。

## Design Philosophy / 设计哲学

> "From Pipeline Encoding to Judgment Encoding"
> "从流程编码到判断编码"

Traditional pipelines encode *what commands to run*. Fan-skill encodes *how to decide what to run*.
传统 pipeline 编码"跑什么命令"。Fan-skill 编码"如何决定跑什么"。

## PR Process / PR 流程

1. Fork the repository / Fork 仓库
2. Create a branch: `git checkout -b add-<analysis-name>` / 创建分支
3. Create 4 files in `knowledge-base/<name>/` / 创建 4 个文件
4. Run `engine/validate_entry.sh` / 运行验证
5. Commit and push / 提交推送
6. Open a Pull Request / 发起 PR
7. CI will automatically validate your entry / CI 自动验证
8. Review and merge / 审核合并
