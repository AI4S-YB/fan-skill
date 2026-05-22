# plant-bioinfo-skills

植物农业组学数据分析 Skill 库 — 面向 Claude Code 的决策引导型生物信息学分析工具集。

## 设计理念

本库不编码固定的分析步骤，而是编码**分析师的判断能力**：
- **Rule 模式**：结构化决策矩阵，条件 → 方法 → 参数
- **Expert 模式**：专家分析笔记，引导 Agent 根据数据特征自主推理
- **Hybrid 模式**（默认）：规则优先，无匹配时自动切换到专家推理

## 包含的 Skill

| Skill | 分析类型 | 状态 |
|-------|---------|:----:|
| bio-plant-infra | 基础设施（数据探查/环境自检/物种速查） | ✅ |
| bio-plant-gwas | 全基因组关联分析 | ✅ |
| bio-plant-population | 群体遗传结构分析 | ✅ |
| bio-plant-genomic-selection | 基因组选择/预测 | 🔜 |
| bio-plant-rnaseq | 转录组差异表达 | 🔜 |
| bio-plant-comparative | 比较基因组学 | 🔜 |
| bio-plant-marker | 分子标记开发 | 🔜 |

## 快速开始

```bash
# 1. 安装
bash install.sh

# 2. 准备输入文件
# samplesheet.csv — 样本→数据映射
# params.yaml — 运行参数 + 决策模式

# 3. 开始分析
# 在 Claude Code 中: "用 bio-plant-gwas 分析我的数据"
```

## 决策模式

在 `params.yaml` 中设置：
```yaml
decision_mode: hybrid  # rule | expert | hybrid (推荐)
```

## 物种支持

水稻、玉米、小麦、大豆、棉花、油菜、拟南芥、番茄、马铃薯、大麦、高粱、甘蔗

## 许可

MIT
