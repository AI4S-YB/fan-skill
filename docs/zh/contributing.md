# 贡献 Fan-Skill

感谢你为 fan-skill 做出贡献。本指南说明项目的结构以及如何添加新的分析能力。

## 架构：B+C 模式

Fan-skill 采用 **B+C（基础 + 配置）** 架构。每个分析条目由两层组成：

- **B（Base，基础层）** — 一个可直接运行的分析笔记本（`.md` 文件），包含完整的 R/Python 代码、植物特有的设计门控、QC 检查点以及故障排除钩子。这是调用 skill 时实际执行的内容。
- **C（Configuration，配置层）** — 一个 `rules.yaml` 文件，描述何时以及如何使用该分析，包括意图匹配、输入验证、设计门控逻辑以及与其他条目的集成。

这种分离意味着领域专家编写分析代码（B），而 skill 的推理引擎使用配置（C）来为用户上下文选择、定制和执行正确的分析。

## 项目结构

```
fan-skill/
├── knowledge-base/          # B+C 条目，每个分析一个目录
│   ├── gwas/                # 示例条目
│   │   ├── rules.yaml       # C: 配置
│   │   └── notebook.md      # B: 分析笔记本
│   ├── rnaseq/
│   ├── genomic-selection/
│   └── ... (27+ 项)
├── tool-registry/           # 工具封装（PLINK、GAPIT、BWA 等）
├── references/             # 参考数据（物种速查表、QC 阈值）
├── templates/              # 可复用文档模板
├── engine/                 # 核心引擎脚本（意图匹配、环境检查）
├── docs/                   # 文档（当前目录）
│   ├── en/                 # 英文文档
│   └── zh/                 # 中文文档
└── theme/                  # 输出样式
```

## 添加新分析：4 个文件

要贡献新的分析能力，请创建以下 4 个文件：

### 1. `rules.yaml` — 配置

```yaml
# knowledge-base/<条目名称>/rules.yaml
entry: <条目名称>
display_name: <人类可读的名称>
category: <能力目录中的某个分类>

intent:
  triggers:             # 触发该条目的关键词/短语
    - <触发短语 1>
    - <触发短语 2>
  anti_triggers:        # 抑制该条目的关键词
    - <反触发词>

inputs:
  required:
    - name: <输入名称>
      type: <文件类型>
      description: <说明>
      format: [<接受的格式>]
  optional:
    - name: <输入名称>
      type: <文件类型>
      description: <说明>

design_gates:           # 植物特有的实验设计检查
  - gate: <门控名称>
    check: <检查内容>
    action: <检查失败时的处理方式>

outputs:
  - <输出描述 1>
  - <输出描述 2>

depends_on:             # 可选：该条目依赖的其他条目
  - <条目名称>

related:                # 可选：相关的其他条目
  - <条目名称>
```

### 2. `notebook.md` — 分析笔记本

笔记本是一个将解释与可执行代码相结合的文学化编程文档。它必须包含：

- **目标**：该分析做什么以及何时使用
- **植物特有的设计门控**：运行前检查倍性、育种系统、群体结构
- **输入清单**：使用格式检查验证所有输入
- **分步分析**：带有内联解释的 R 或 Python 代码块
- **QC 检查点**：在关键步骤验证中间结果（如 MAF 分布、PCA 异常值）
- **输出表格和图表**：发表级别的可视化图表和格式化表格
- **故障排除钩子**：常见失败模式及恢复方法

代码块使用标准的围栏格式：

````markdown
```r
# R 代码
library(ggplot2)
...
```
````

### 3. `tool-registry/<工具名称>.md` — 工具封装（如需要）

如果你的分析使用的工具尚未在工具注册表中，请创建一个封装：

```markdown
# <工具名称>

## 安装
...

## 参数参考
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|

## 参数决策表
| 场景 | 推荐设置 | 理由 |
|------|----------|------|

## 使用
...
```

### 4. 更新 `capability-catalog.md`

将你的条目添加到 `docs/en/capability-catalog.md` 和 `docs/zh/capability-catalog.md` 中相应的分类表格。

## 设计门控模式

设计门控是 fan-skill 的一个关键特性。它们编码了植物特有的实验设计知识，防止执行不合适的分析。每个笔记本必须在开头实现设计门控。

示例：在对作物物种运行 GWAS 之前，检查：

1. **倍性门控**：是否为同源多倍体？如果是，将基因型重新编码为二倍体或使用多倍体感知方法。
2. **群体结构门控**：是否存在未控制的群体结构？如果是，纳入 PCA 协变量。
3. **自交/异交门控**：自交作物 → 使用 FarmCPU/BLINK。异交作物 → 考虑使用带亲缘关系的 MLM。
4. **LD 衰减门控**：估计 LD 衰减距离以调整显著性阈值。

每个门控在 `rules.yaml` 中有三个组成部分：

```yaml
design_gates:
  - gate: <门控名称>
    check: <检查内容>
    action: <检查失败时的处理方式>
```

以及在 `notebook.md` 中的对应实现：

```r
# ---- 设计门控：群体结构 ----
# 检查：是否存在未控制的群体结构？
pca <- read.table("pca.eigenvec", header=FALSE)
var_explained <- pca_eigenval / sum(pca_eigenval) * 100
if (var_explained[1] > 10) {
  message("警告：PC1 解释 >10% 方差。")
  message("操作：纳入前几个 PC 作为协变量。")
  covariates <- pca[, 1:3]  # 使用前 3 个 PC
}
```

## Pull Request 流程

1. **Fork** 仓库并创建功能分支
2. **创建上述 4 个文件**，对应你的新分析
3. **验证** 你的 `rules.yaml` 是否符合模式：
   ```bash
   python engine/validate_rules.py knowledge-base/<条目名称>/rules.yaml
   ```
4. **测试** 你的笔记本，使用 `test_matrix.csv` 中提供的示例数据进行端到端运行
5. **检查集成问题**：你的条目不应与现有条目的意图触发词冲突
6. **提交 PR**，描述内容包括：
   - 你添加的是什么分析
   - 你实现了哪些植物特有的设计门控
   - 任何新的工具注册表条目
   - 测试结果

## 风格指南

- **笔记本代码**：使用 `set.seed()` 以保证可重现性。R 代码注释使用英文，markdown 叙述可用中文或英文。
- **rules.yaml**：所有键和值使用英文（推理引擎通过程序解析该文件）。
- **设计门控**：门控应为建设性的（建议修复方案），而非阻塞性的（无替代方案直接拒绝）。
- **参考文献**：在适当的地方引用方法论文。使用标准期刊缩写。

## 有问题？

在 GitHub 上提 issue 或发起讨论。我们特别欢迎针对非模式植物物种、多倍体特有方法以及与公共植物数据库（Ensembl Plants、Phytozome、Gramene 等）集成的贡献。
