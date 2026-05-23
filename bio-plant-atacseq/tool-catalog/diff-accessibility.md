# 差异可及性分析 -- 工具目录

## 概述

差异可及性分析用于鉴定不同条件之间（如处理 vs 对照、不同发育阶段）开放染色质状态的显著变化区域。

## 推荐工具

### 1. DiffBind + DESeq2 -- 强烈推荐

**描述**: DiffBind 是差异结合/可及性分析的标准工具，使用 DESeq2 或 edgeR 作为统计引擎。需要至少 2 个生物学重复。

**R 环境设置**:
```r
library(DiffBind)
library(DESeq2)
```

**完整分析流程**:
```r
# 步骤 1: 准备 samplesheet
# samplesheet.csv 格式:
# SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller
# sample1,leaf,ATAC,control,1,data/sample1.bam,data/sample1_peaks.narrowPeak,narrow
# sample2,leaf,ATAC,control,2,data/sample2.bam,data/sample2_peaks.narrowPeak,narrow
# sample3,leaf,ATAC,treatment,1,data/sample3.bam,data/sample3_peaks.narrowPeak,narrow
# sample4,leaf,ATAC,treatment,2,data/sample4.bam,data/sample4_peaks.narrowPeak,narrow

samples <- read.csv("samplesheet.csv")
dba_obj <- dba(sampleSheet = samples)

# 步骤 2: 计算 overlap
dba_obj <- dba.peakset(dba_obj,
  consensus = DBA_CONDITION,
  minOverlap = 2)
```

**参数说明**:
- `DBA_CONDITION`: 差异分析的分组依据
- `DBA_DESEQ2`: 使用 DESeq2 统计方法（推荐）
- `minMembers = 2`: 每个分组至少 2 个样本
- `bUseSummarizeOverlaps = TRUE`: 使用更准确的计数方法

```r
# 步骤 3: 计数
dba_obj <- dba.count(dba_obj,
  bUseSummarizeOverlaps = TRUE,
  score = DBA_SCORE_READS)

# 步骤 4: 差异分析
dba_obj <- dba.contrast(dba_obj,
  categories = DBA_CONDITION,
  minMembers = 2)
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

# 步骤 5: 提取结果
results <- dba.report(dba_obj,
  method = DBA_DESEQ2,
  th = 1,         # FDR threshold
  fold = 0,        # no fold-change cutoff
  bCounts = TRUE)  # include read counts

# 导出
write.csv(as.data.frame(results),
  "differential_accessibility.csv")
```

### 2. 归一化方法选择

对 ATAC-seq 数据，推荐以下归一化方法：

**建议（按顺序尝试）**:
1. **DESeq2 默认**（median-of-ratios）: 适用于大多数情况
2. **RLE (Relative Log Expression)**: 适用于测序深度差异大的情况
3. **TMM (Trimmed Mean of M-values)**: 备选方法

```r
dba_obj <- dba.normalize(dba_obj,
  method = DBA_DESEQ2,
  normalize = DBA_NORM_DEFAULT)
```

### 3. 差异可及性阈值

```r
# 显著差异可及区域 (DARs)
results_sig <- as.data.frame(results)
results_sig <- results_sig[results_sig$FDR < 0.05, ]

# 增加变化幅度 (至少 2 倍变化)
results_sig_fc <- results_sig[abs(results_sig$Fold) > 1, ]
```

**统计总结**:
```r
# 统计 DAR 数目
dba.show(dba_obj, bContrasts = TRUE)

cat("Total significant DARs (FDR < 0.05):", nrow(results_sig), "\n")
cat("FDR < 0.05 & |FC| > 1:", nrow(results_sig_fc), "\n")
```

---

## 可视化

### MA Plot
```r
dba.plotMA(dba_obj, contrast = 1)
```

### PCA Plot
```r
dba.plotPCA(dba_obj,
  attributes = c(DBA_CONDITION, DBA_REPLICATE),
  vColors = c("blue", "red"))
```

### Volcano Plot
```r
dba.plotVolcano(dba_obj)
```

### 相关性热图
```r
dba.plotHeatmap(dba_obj, contrast = 1,
  correlations = TRUE, colScheme = "Reds")
```

---

## 注意事项

1. **生物学重复**: 必须至少 2 个重复才能使用 DESeq2
2. **Consensus peaks**: 在计数之前，先确定每个条件下至少 2 个样本共有的 consensus peaks
3. **批次效应**: 如果 PCA 显示批次效应而非生物学差异，考虑使用 Combat-seq 或 RUVseq
4. **低重复统计**: 如果只有 2 个重复，DESeq2 可能过于保守；可考虑使用较宽松的阈值
5. **多倍体植物**: 使用 featureCounts 的 multi-mapping 选项处理同源区域
