# 常见分析陷阱

**注意**: 分析特定的陷阱和避免策略现已分散到各知识条目的 `notebook.md` 中。

## 如何使用

每个分析类型的 B 层专家推理笔记 (`notebook.md`) 中包含：
- 该分析类型专属的常见陷阱 (Common Pitfalls)
- 陷阱产生的原因和机制
- 如何检测和避免的建议

当需要了解某个特定分析的陷阱时，请查阅对应条目的 `notebook.md`。

## 陷阱位置

| 分析类型 | 陷阱在 notebook.md 中 |
|---------|---------------------|
| GWAS | `knowledge-base/gwas/notebook.md` |
| RNA-seq | `knowledge-base/rnaseq/notebook.md` |
| 基因组选择 | `knowledge-base/genomic-selection/notebook.md` |
| 比较基因组学 | `knowledge-base/comparative/notebook.md` |
| 扩增子测序 | `knowledge-base/amplicon/notebook.md` |
| ATAC-seq | `knowledge-base/atacseq/notebook.md` |
| ChIP-seq | `knowledge-base/chipseq/notebook.md` |
| CRISPR | `knowledge-base/crispr/notebook.md` |
| Enviromics | `knowledge-base/enviromics/notebook.md` |
| eQTL | `knowledge-base/eqtl/notebook.md` |
| 基因组注释 | `knowledge-base/genome-annotation/notebook.md` |
| 基因组组装 | `knowledge-base/genome-assembly/notebook.md` |
| 基因型填充 | `knowledge-base/genotype-imputation/notebook.md` |
| GRN | `knowledge-base/grn/notebook.md` |
| 杂种预测 | `knowledge-base/hybrid-prediction/notebook.md` |
| 分子标记 | `knowledge-base/marker/notebook.md` |
| 代谢组学 | `knowledge-base/metabolomics/notebook.md` |
| 宏基因组学 | `knowledge-base/metagenomics/notebook.md` |
| 甲基化组学 | `knowledge-base/methylation/notebook.md` |
| 多组学整合 | `knowledge-base/multi-omics/notebook.md` |
| 泛基因组 | `knowledge-base/pan-genome/notebook.md` |
| 表型分析 | `knowledge-base/phenotype/notebook.md` |
| 群体遗传 | `knowledge-base/population/notebook.md` |
| 蛋白质组学 | `knowledge-base/proteomics/notebook.md` |
| QTL 作图 | `knowledge-base/qtl-mapping/notebook.md` |
| 小 RNA | `knowledge-base/small-rna/notebook.md` |
| 时间序列 | `knowledge-base/time-series/notebook.md` |
| 变异检测 | `knowledge-base/variant-calling/notebook.md` |
| 可视化 | `knowledge-base/visualization/notebook.md` |

## 通用陷阱（跨分析类型）

以下陷阱适用于几乎所有分析类型，保留在此以供快速参考：

### 生物学重复不足
任何涉及统计推断的分析都需要生物学重复。2 个重复 → 方差估计极不稳定，统计功效很低。3 个重复是共识最低标准。

### 缺对照组
没有对照组，无法区分"处理效应"和"本来就有（或技术噪声引起）"的差异。

### 批次效应
不同批次（建库日期、操作人员、测序 lane）之间可能产生系统性偏差，掩盖真实的生物学信号。始终在分析前检查 PCA/热图看样本是否按批次而非生物学分组聚类。

### 多重检验校正
全基因组尺度的检验（GWAS 数百万 SNP，RNA-seq 数万基因）必须进行多重检验校正。未校正的 p < 0.05 会产生大量假阳性。

**新增分析类型时，在该条目的 `notebook.md` 中记录其专属陷阱即可，无需修改此文件。**
