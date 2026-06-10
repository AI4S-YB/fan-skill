# fan-skill

AI 驱动的植物生物信息与育种分析引擎。从生物学问题到发表级结果。三层知识架构（个人 → 通用 → Agent 兜底）配合 B+C 双模式决策系统。

> 🇬🇧 [English Docs](README.md)

## 快速开始

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

在 Claude Code 中直接说出你的生物学问题：

> "我有 300 份水稻材料的 GBS 数据和 3 年的产量表型，想找到控制粒重的基因位点"

Fan-skill 会自动匹配知识库、检查实验设计、推荐分析路径并执行分析。

📖 **[快速开始指南](docs/zh/quick-start.md)** · **[分析能力目录](docs/zh/capability-catalog.md)** · **[贡献指南](docs/zh/contributing.md)**

## 架构

```
用户说出生物学问题
        │
        ▼
┌─────────────────────────────────┐
│  SKILL.md                       │  ← 统一入口
│  意图识别 → 三层匹配 → 执行
└──────────────┬──────────────────┘
               │
    ┌──────────┼──────────┐
    ▼          ▼          ▼
┌────────┐ ┌───────┐ ┌───────┐
│  用户  │ │ gwas  │ │rnaseq │  ← user-knowledge/ (你的经验)
│  条目  │ │ rules │ │ rules │     knowledge-base/ (31 条目)
└───┬────┘ └───┬───┘ └───┬───┘
    │          │         │
    └──────────┼─────────┘
               ▼
┌─────────────────────────────────┐
│  engine/                        │  ← 共享执行引擎
│  rule_engine · validate · log   │
└─────────────────────────────────┘
               │
    ┌──────────┼──────────┐
    ▼          ▼          ▼
  分析报告   决策日志   图表
               │
    (未匹配时) Agent 兜底推理 → 沉淀到 user-knowledge/
```

经过苹果 RNA-seq 和水稻候选基因分析真实验证。每项决策可追溯。

### 三层知识架构

| 层次 | 路径 | 优先级 | 内容 |
|------|------|:------:|------|
| **用户** | `user-knowledge/` | 🔺 最高 | 你的个人分析经验。初始为空，越用越丰富。 |
| **通用** | `knowledge-base/` | 🔸 默认 | 31 条精选条目，覆盖主流植物生信分析。 |
| **Agent 兜底** | (模型推理) | 🔻 补充 | 两层均未匹配时，Agent 自主推理并可将结果沉淀到用户层。 |

### B+C 双模式

| 层次 | 作用 | 形式 |
|------|------|------|
| **C 层 (Rule)** | 确定性决策 — 同样数据永远选同样方法 | `rules.yaml` |
| **B 层 (Expert)** | 灵活推理 — 规则未覆盖时像专家一样判断 | `notebook.md` |
| **用户控制** | `decision_mode: rule \| expert \| hybrid` | `params.yaml` |
| **用户条目** | 最小结构 `meta.yaml` + `notebook.md`，可选 `rules.yaml` | `user-knowledge/` |

## 分析能力

| 目标 | fan-skill 可以 |
|------|-------------|
| 找控制性状的基因 | GWAS · QTL 定位 · eQTL · 候选基因关联 |
| 理解基因功能 | RNA-seq 差异表达 · 时间序列 · GRN 调控网络 · 多组学 |
| 预测育种值 | 基因组选择 · 杂种优势预测 · 环境组学 |
| 数据预处理 | 变异检测 · 基因型填充 · 表型分析 |
| 表观组学 | ChIP-seq · ATAC-seq · DNA 甲基化 |
| 微生物组 | 16S/ITS 扩增子 · 宏基因组 |
| 代谢与蛋白组 | LC-MS/GC-MS · DIA/DDA · PPI |
| 进化分析 | 比较基因组 · 泛基因组 |
| 编辑与育种 | CRISPR 设计 · 分子标记开发 |
| 科研图表 | 发表级 ggplot2 图表 |

完整目录：**[分析能力目录](docs/zh/capability-catalog.md)**

## 设计理念

Fan-skill 不编码分析步骤，而是编码**分析师的判断能力**——即专家如何决定"该做什么分析、该用什么方法"。

- **Pipeline 的问题**：固定步骤无法适应多样化的数据
- **原子 Skill 的问题**：自由组合导致调用路径不稳定
- **Fan-skill 的方案**：三层知识架构。你的个人经验（`user-knowledge/`）优先于 31 条精选条目（`knowledge-base/`）。当两层都不匹配时，Agent 推理作为兜底，并可保存为新的个人条目——系统会越用越强。
- **B+C 双模式**：C 层编码"选什么方法"的结构化规则，B 层编码"为什么"的专家推理。两种模式在每个条目中共存。

每项决策可追溯——`engine/log_decision.sh` 记录每一步的选择依据。

## 贡献

**贡献到通用知识库：** 新增分析只需 4 个文件（不是 14 个）。详见 **[贡献指南](docs/zh/contributing.md)**。

```bash
cp templates/rules-template.yaml knowledge-base/<名称>/rules.yaml
cp templates/notebook-template.md knowledge-base/<名称>/notebook.md
# + consult-guide.md + analysis-primer.md
engine/validate_entry.sh knowledge-base/<名称>/
```

**积累个人知识库：** 正常使用即可。当 Agent 兜底解决了一个新分析时，选择保存——条目自动生成在 `user-knowledge/drafts/`。审核确认后移到 `user-knowledge/confirmed/`。

## 许可

MIT
