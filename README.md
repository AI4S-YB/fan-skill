# fan-skill

AI 驱动的植物生物信息与育种分析引擎。基于知识库的 B+C 双模式决策架构，从生物学问题到发表级结果。

## 快速开始

```bash
git clone git@github.com:AI4S-YB/fan-skill.git
cd fan-skill
bash install.sh
```

在 Claude Code 中直接说出你的问题：

> "我有 300 份水稻的 GBS 数据和 3 年的产量表型，想找到控制产量的基因位点"

fan-skill 会引导对话、匹配知识库、推荐分析路径、执行分析、生成报告。

详见：**[Quick Start / 快速开始](docs/quick-start.md)** · **[Analysis Capabilities / 分析能力](docs/capability-catalog.md)** · **[Contributing / 贡献指南](docs/contributing.md)**

## 架构

```
用户说出生物学问题
        │
        ▼
┌─────────────────────────────────┐
│  SKILL.md                       │  ← 统一入口 (大脑)
│  意图识别 → 知识匹配 → 路径发现  │
└──────────────┬──────────────────┘
               │
    ┌──────────┼──────────┐
    ▼          ▼          ▼
┌───────┐ ┌───────┐ ┌───────┐
│ gwas  │ │rnaseq │ │  ...  │  ← knowledge-base/ (领域知识)
│ rules │ │ rules │ │       │    每个条目: C层决策规则 + B层专家笔记
└───┬───┘ └───┬───┘ └───┬───┘
    │         │         │
    └─────────┼─────────┘
              ▼
┌─────────────────────────────────┐
│  engine/                        │  ← 共享执行引擎
│  rule_engine · discover_chains  │
│  run_pipeline · log_decision    │
└─────────────────────────────────┘
              │
              ▼
        分析报告 + 决策日志 + 图表
```

### B+C 双模式

| 层次 | 作用 | 形式 |
|------|------|------|
| **C 层 (Rule)** | 确定性决策 — 同样的数据特征永远选同样的方法 | `rules.yaml` |
| **B 层 (Expert)** | 灵活推理 — 规则未覆盖时像专家一样判断 | `notebook.md` |
| **用户控制** | `decision_mode: rule \| expert \| hybrid` | `params.yaml` |

## 能力

### 分析场景

| 你想做什么 | fan-skill 可以 |
|-----------|---------------|
| 找控制性状的基因 | GWAS · QTL 定位 · eQTL · 候选基因注释 |
| 预测育种值、选配亲本 | 基因组选择 · 配合力预测 · 最优杂交组合 |
| 理解基因功能 | RNA-seq 差异表达 · 时间序列 · GRN 调控网络 |
| 看群体遗传结构 | PCA · ADMIXTURE · Fst · 系统发育树 |
| 进化分析 | 共线性 · Ks 分析 · 基因家族 · 选择压力 |
| 多组学整合 | 转录+代谢+蛋白联合 · 多组学因子分析 |
| 开发检测标记 | KASP · InDel · SSR · 功能性标记 · 亲本推荐 |
| 表观组学 | ChIP-seq · ATAC-seq · DNA 甲基化 |
| 基因组组装与注释 | 长读长组装 · 重复序列屏蔽 · 基因预测 |
| 微生物组 | 16S/ITS 扩增子 · 宏基因组 binning |
| 代谢组与蛋白组 | LC-MS/GC-MS · DIA/DDA · PPI 网络 |
| CRISPR 设计 | sgRNA 设计 · 脱靶预测 |
| 数据预处理 | 变异检测 · 基因型填充 · 表型 BLUP/遗传力 |
| 科研图表 | 发表级 ggplot2 · 色盲友好 · 矢量输出 |

完整能力列表见 `knowledge-base/` 目录。

## 设计理念

fan-skill 不编码分析步骤，而是编码**分析师的判断能力**：

- **Pipeline 模式的问题**：固定的分析步骤无法适应多样化的数据。"2 个重复的 RNA-seq"和"300 份材料的 GWAS"需要完全不同的处理。
- **原子 Skill 的问题**：自由组合导致调用路径不稳定，同样的任务在不同 session 中可能走不同的路径。
- **fan-skill 的方案**：知识条目的 C 层编码"该用什么方法"的结构化规则，B 层编码"为什么"的专家推理。用户通过 `decision_mode` 控制稳定性与灵活性的平衡。

每项分析决策可追溯——`engine/log_decision.sh` 记录每一步的选择依据。

## 贡献

新增分析能力只需创建 **2 个文件**（而非传统的完整 Skill 目录）：

```
knowledge-base/<your-analysis>/
├── rules.yaml      ← C 层：决策规则 + inputs/outputs 声明
└── notebook.md     ← B 层：专家分析笔记
```

1. 从 `templates/` 复制模板
2. 填写领域知识
3. 运行 `engine/validate_entry.sh knowledge-base/<your-analysis>/`
4. 提交 PR

详见 `templates/CONTRIBUTING.md`。

## 物种支持

水稻、玉米、小麦、大豆、棉花、油菜、拟南芥、番茄、马铃薯、大麦、高粱、甘蔗

## 许可

MIT
