# ATAC-seq 可视化 -- 工具目录

## 概述

ATAC-seq 可视化覆盖多个层面：基因组浏览器 tracks、peak 分布概况、TSS 可及性谱、差异可及性热图和 TF footprint 可视化。

## 推荐工具

### 1. deepTools -- 核心可视化工具

**描述**: deepTools 是 ChIP/ATAC-seq 数据可视化的标准工具集，支持连续信号的生成和多种分析图。

**bigWig 生成**:
```bash
# CPM 归一化的 bigWig track
bamCoverage -b sample_dedup.bam \
  --normalizeUsing CPM \
  --binSize 10 \
  --effectiveGenomeSize 119481543 \
  -o sample_CPM.bw \
  -p 8

# RPKM 归一化
bamCoverage -b sample_dedup.bam \
  --normalizeUsing RPKM \
  --binSize 10 \
  -o sample_RPKM.bw \
  -p 8
```

**TSS 可及性谱**:
```bash
# TSS 区域信号矩阵
computeMatrix reference-point \
  -S sample1_CPM.bw sample2_CPM.bw \
  -R genes.bed \
  --referencePoint TSS \
  -a 2000 -b 2000 \
  --binSize 50 \
  -o matrix_TSS.gz \
  -p 8

# TSS 信号分布图
plotProfile -m matrix_TSS.gz \
  -o profile_TSS.pdf \
  --perGroup \
  --plotTitle "ATAC-seq TSS Enrichment" \
  --yAxisLabel "CPM"

# TSS 热图
plotHeatmap -m matrix_TSS.gz \
  -o heatmap_TSS.pdf \
  --colorMap YlOrRd \
  --zMax 10 \
  --whatToShow "plot, heatmap and colorbar" \
  --yAxisLabel "Genes" \
  --xAxisLabel "Distance from TSS (bp)"
```

### 2. IGV (Integrative Genomics Viewer) -- 基因组浏览器

**描述**: 交互式可视化基因组上的 ATAC-seq 信号。

加载内容:
- `sample_CPM.bw` -- ATAC-seq 连续信号
- `sample_peaks.narrowPeak` -- called peaks 位置
- `genome.gff` -- 基因注释

### 3. 插入片段长度分布

```r
library(ggplot2)

# 读取片段长度
frag_sizes <- read.table("fragment_sizes.txt",
  col.names = c("count", "length"))

# 绘图
ggplot(frag_sizes, aes(x = length, y = count)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(147, 294, 441),
    linetype = "dashed", color = "red") +
  xlim(0, 600) +
  labs(title = "ATAC-seq Fragment Length Distribution",
       x = "Fragment Length (bp)",
       y = "Count") +
  theme_minimal()
```

**解读**:
- **< 150bp**: Nucleosome-free region (NFR)，来自开放染色质
- **180-247bp**: Mono-nucleosome，一个核小体保护
- **315-473bp**: Di-nucleosome，两个核小体
- **440-680bp**: Tri-nucleosome

### 4. 可及性峰注释分布饼图

```r
library(ChIPseeker)

# 生成基因组注释图
peak_anno <- annotatePeak(peaks, TxDb = txdb,
  tssRegion = c(-3000, 3000))

# 饼图
plotAnnoPie(peak_anno,
  main = "ATAC-seq Peak Genomic Distribution")

# 到 TSS 距离分布
plotDistToTSS(peak_anno,
  title = "Peak Distribution Relative to TSS")
```

### 5. 差异可及性可视化

```r
library(DiffBind)

# MA plot
dba.plotMA(dba_obj, contrast = 1,
  main = "MA Plot: Treatment vs Control")

# PCA
dba.plotPCA(dba_obj,
  attributes = c(DBA_CONDITION, DBA_REPLICATE),
  label = DBA_ID)

# 火山图
dba.plotVolcano(dba_obj, contrast = 1,
  main = "Volcano: Treatment vs Control")

# 相关性热图
dba.plotHeatmap(dba_obj, contrast = 1,
  correlations = TRUE,
  colScheme = "Reds")
```

### 6. TF Footprint 可视化

```r
# TOBIAS plot 函数
library(ggplot2)

# Footprint profile
ggplot(footprint_data, aes(x = position, y = signal)) +
  geom_line(aes(color = condition), size = 1) +
  geom_vline(xintercept = c(0, motif_length),
    linetype = "dashed") +
  labs(x = "Position relative to motif",
       y = "Tn5 insertion signal",
       title = "TF Footprint Profile") +
  theme_minimal()
```

### 7. Heatmap + Metagene Plots (EnrichedHeatmap)

```r
library(EnrichedHeatmap)
library(circlize)

# 信号矩阵
mat <- normalizeToMatrix(sample_CPM.bw, peaks,
  value_column = "score",
  extend = 2000, mean_mode = "w0",
  w = 50)

# 热图
EnrichedHeatmap(mat,
  col = colorRamp2(c(0, 5, 10),
    c("white", "orange", "red")),
  name = "CPM",
  top_annotation = HeatmapAnnotation(
    enriched = anno_enriched(
      ylim = c(0, 10))))
```

---

## 植物特有可视化考虑

1. **染色体命名**: 植物染色体命名不一致，需在 plots 中使用规范名称
2. **大型基因组**: 玉米/小麦等大型基因组在基因组浏览器中可能加载缓慢，使用 100bp bin 而非 10bp bin
3. **多倍体基因组**: 亚基因组特异的染色体应分开标注，如 A01, B01, D01 格式
4. **叶绿体噪声**: 如果是绿色组织 ATAC-seq，可选择性展示核基因组区域
5. **植物特殊区域可视化**: 着丝粒区/端粒区可能显示异常高信号（偏倚），应做适当注释

---

## 质量报告生成

```bash
# MultiQC 汇总所有 QC 指标
multiqc qc/raw/ qc/trimmed/ star_logs/ \
  -o multiqc_report/ \
  -n atacseq_multiqc_report
```
