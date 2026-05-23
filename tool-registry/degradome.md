# 降解组分析 — 工具目录

## 概述

降解组测序 (Degradome-seq / PARE) 用于验证 miRNA 对靶基因的剪切。通过检测 miRNA 引导的 mRNA 剪切产物，可以高通量地确认 miRNA 靶基因。

## 推荐工具

### 1. CleaveLand (推荐)

**描述**: 降解组数据分析的标准工具，用于鉴定 miRNA 引导的靶基因剪切位点。

**工作流程**:
1. 将降解组 reads 比对到转录组
2. 识别 miRNA 靶位点
3. 验证剪切信号 (检测降解组 reads 5‘ 端在 miRNA 第 10-11 位碱基对应的位置富集)

**输入**:
- 降解组测序 reads (FASTA)
- miRNA 序列 (FASTA)
- 转录组序列 (FASTA)

**输出**:
- 剪切位点分类 (Category 0-4)
- 靶基因 T-plot
- 验证的 miRNA-靶基因对列表

**参考**: https://github.com/MikeAxtell/CleaveLand4

### 2. PAREsnip2 (备选)

**描述**: 基于网页的降解组分析工具，提供交互式可视化。

**适用场景**: 快速交互式分析，无需命令行操作。
