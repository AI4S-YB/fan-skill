# 差异结合分析 — 工具目录

## 概述

差异结合分析 (Differential Binding Analysis) 用于识别不同实验条件之间蛋白质-DNA 结合强度显著变化的基因组区域。

## 推荐工具

### 1. DiffBind (推荐)

**描述**: 用于 ChIP-seq 差异结合分析的 R/Bioconductor 包，集成了 DESeq2 和 edgeR 统计方法。

**工作流程**:
1. 读取所有样本的 peak 集合
2. 合并 peaks 为一致的 peak set (consensus peaks)
3. 计算每个 peak 区域的 reads 计数
4. 使用 DESeq2 或 edgeR 进行差异分析
5. 结果可视化和导出

**R 代码示例**:
```r
library(DiffBind)

# 读取样本信息表
samples <- read.csv("samplesheet.csv")

# 创建 DBA 对象
dba_obj <- dba(sampleSheet = samples)

# 计算计数矩阵
dba_obj <- dba.count(dba_obj, bUseSummarizeOverlaps = TRUE)

# 设置对比
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION,
                         minMembers = 2)

# 差异分析
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

# 提取结果
results <- dba.report(dba_obj, method = DBA_DESEQ2)
```

**样本信息表格式 (samplesheet.csv)**:
| SampleID | Tissue | Factor | Condition | Replicate | bamReads       | Peaks              | PeakCaller |
|----------|--------|--------|-----------|-----------|----------------|--------------------|------------|
| WT_1     | leaf   | TF1    | WT        | 1         | wt1_sorted.bam | wt1_peaks.narrowPeak | narrow     |
| WT_2     | leaf   | TF1    | WT        | 2         | wt2_sorted.bam | wt2_peaks.narrowPeak | narrow     |
| Mut_1    | leaf   | TF1    | mutant    | 1         | mut1_sorted.bam| mut1_peaks.narrowPeak| narrow     |
| Mut_2    | leaf   | TF1    | mutant    | 2         | mut2_sorted.bam| mut2_peaks.narrowPeak| narrow     |

**输出**:
- 差异结合的基因组区域列表 (BED/GFF)
- MA plot, Volcano plot
- PCA plot

### 2. DESeq2 (直接使用)

**描述**: 如果已有 peak 区域的 count matrix，可以直接使用 DESeq2 进行分析。

**关键参数**:
- `design = ~ condition`: 实验设计
- `contrast = c("condition", "treatment", "control")`: 对比设置

### 3. edgeR (备选)

**描述**: 适用于样本数较少的情况。

---

## 植物特殊注意事项

1. **样本复杂度**: 植物组织中的细胞壁和次生代谢物可能影响 ChIP 效率，需要足够 Input 对照
2. **重复数要求**: 推荐至少 2 个生物学重复，3 个更好
3. **peak 一致性**: 植物 ChIP-seq 数据噪音可能较大，注意 consensus peak 的选取标准
4. **多倍体处理**: 同源基因区域的 reads 分配可能不准确，影响差异分析
