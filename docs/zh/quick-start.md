# Fan-Skill 快速开始

## 安装

```bash
git clone git@github.com:AI4S-YB/fan-skill.git
cd fan-skill

# 选择你的 AI 编程 Agent 安装:
bash install-claude.sh     # Claude Code
bash install-codex.sh      # Codex CLI
bash install-gemini.sh     # Gemini CLI
bash install-opencode.sh   # OpenCode
bash install-all.sh        # 自动检测并全部安装

# 或使用通用安装器:
# npx skills add AI4S-YB/fan-skill
```

## 第一次使用

在 Claude Code 中直接用自然语言描述你的生物学问题：

> "我有 300 份水稻材料的 GBS 数据和 3 年的产量表型，想找到控制粒重的基因位点"

Fan-skill 会自动匹配知识库、检查实验设计、推荐分析路径并执行分析。

## 使用方式

| 模式 | 场景 | 示例 |
|------|------|------|
| **我知道要做什么** | 有数据，有明确的分析目标 | "对我的 VCF 和表型数据跑 GWAS" |
| **我有数据，能做什么？** | 有数据但不确定能做什么分析 | "我有干旱处理水稻的 RNA-seq 数据，能分析什么？" |
| **我有生物学问题** | 想回答研究问题但不知道用什么分析 | "找到控制耐盐的基因" |

## 环境要求

- Claude Code（任何较新版本）
- R 4.0+（大部分分析需要）
- Python 3.10+（数据探查和 GFF3 解析）
- 个别分析可能需要额外工具（PLINK、GAPIT、BWA 等）。运行 `engine/check_env.sh` 检查。

## 下一步

- 查看全部 30 项分析能力：[分析能力目录](capability-catalog.md)
- 贡献新分析：[贡献指南](contributing.md)
- 了解架构：[设计文档](../superpowers/specs/2026-05-22-fan-skill-design.md)
