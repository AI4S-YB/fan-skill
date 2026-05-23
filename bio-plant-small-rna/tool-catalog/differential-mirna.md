# 差异 miRNA 表达分析 — 工具目录

## 概述

差异 miRNA 表达分析旨在识别在不同实验条件 (处理 vs 对照、不同发育阶段等) 下表达水平显著变化的 miRNA。

## 推荐工具

### 1. DESeq2 (推荐)

**描述**: 基于负二项分布的差异表达分析工具，适用于 miRNA 计数数据。

**适用条件**: 至少 3 个生物学重复 (推荐)。

**工作流程**:
1. 准备 miRNA count matrix (raw counts)
2. 创建实验设计矩阵 (treatment vs control)
3. 运行 DESeq2 差异分析
4. 结果筛选 (|log2FC| > 1, padj < 0.05)

**R 包**: `DESeq2`

```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treatment", "control"))
```

**输出**:
- 差异表达 miRNA 列表
- MA plot, Volcano plot
- PCA plot 用于样本质量控制

### 2. edgeR (备选)

**描述**: 适用于样本数较少的情况 (2 vs 2)。

**适用条件**: 2 个生物学重复时可考虑使用。

**R 包**: `edgeR`
