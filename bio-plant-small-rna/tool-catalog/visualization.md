# 小 RNA 可视化 — 工具目录

## 概述

小 RNA 分析结果的可视化对于数据解读和结果展示至关重要。本目录涵盖常用的可视化类型和推荐工具。

## 可视化类型与工具

### 1. 差异表达可视化

**工具**: ggplot2 (R)

**图类型**:
- **Volcano plot**: 展示差异 miRNA 的显著性和变化幅度
- **MA plot**: 展示 log2 fold-change 与平均表达量的关系
- **Heatmap**: 差异 miRNA 的表达模式聚类

**R 代码示例 (Volcano plot)**:
```r
library(ggplot2)
library(ggrepel)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw()
```

### 2. 小 RNA 长度分布

**工具**: ggplot2, Python matplotlib

**描述**: 展示测序 reads 的长度分布 (18-30nt)，反映小 RNA 群体的组成。植物小 RNA 通常以 21nt 和 24nt 为主。

### 3. miRNA 前体二级结构

**工具**: miRDeep2 内置, RNAfold

**描述**: 显示预测的 miRNA 前体二级发夹结构，用于验证 miRNA 预测结果。

### 4. 靶基因 T-plot (降解组分析)

**工具**: CleaveLand 内置

**描述**: 展示 miRNA 引导的靶基因剪切位点，降解组 reads 在剪切位点附近的分布。

### 5. miRNA 基因组分布

**工具**: Circos, IGV

**描述**: 在基因组浏览器中展示 miRNA 的位置、reads 覆盖度等信息。
