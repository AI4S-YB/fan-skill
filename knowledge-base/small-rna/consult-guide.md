# 小 RNA 测序分析引导

## 这个分析解决什么问题

用户想鉴定和分析小 RNA（miRNA、siRNA 等）及其靶基因。典型表述：
- "我想鉴定样品中的 miRNA"
- "哪些 miRNA 在处理条件下差异表达"
- "这个 miRNA 的靶基因是什么"

## 你需要问什么

1. **目标分子**: miRNA / siRNA / piRNA？植物中主要是 miRNA 和 siRNA
2. **物种和数据库**: miRBase 有收录吗？非模式物种需 de novo 预测
3. **实验设计**: 处理/对照？每组几个重复？
4. **分析目标**: 已知 miRNA 鉴定 / 新 miRNA 预测 / 差异表达 / 靶基因预测？
5. **数据质量**: reads 长度分布应在 18-30nt

## 相关信息收集后

收集够信息后 → 进入方法选择 (`rules.yaml`)

## 植物特定考量

- **miRNA 长度分布**: 植物 miRNA 主要分布在 21nt 和 24nt。21nt miRNA 介导 mRNA 切割，24nt miRNA 介导 DNA 甲基化。reads 长度分布可以区分这两类
- **siRNA**: 植物有大量内源 siRNA（特别是 24nt heterochromatic siRNA），通常比 miRNA 丰度更高。需要区分 miRNA 和 siRNA
- **非模式物种**: miRBase 中非模式植物 miRNA 收录有限。需要 de novo 预测（miRDeep2、ShortStack）而非仅依赖已知 miRNA 比对
- **miRNA 靶标预测**: 植物 miRNA-靶基因配对要求高序列互补性（通常 ≤ 4 个错配）。psRNATarget 或 TargetFinder 是常用工具
- **降解组测序**: 如果需要验证 miRNA 是否真的切割靶基因，需要降解组测序（Degradome-seq / PARE）或 5' RACE

## 对话注意事项

- 使用生物学语言而非工具名称解释推荐理由（如"这个分析找控制性状的基因组区域"而非"用 GAPIT 跑混合线性模型"）
- 一次只问一个问题，根据用户的回答自然推进对话
- 如果用户已提供的信息覆盖了上述某个问题，跳过不要重问
- 确保在推荐分析前收集到足够信息——信息不足时推荐不准确
