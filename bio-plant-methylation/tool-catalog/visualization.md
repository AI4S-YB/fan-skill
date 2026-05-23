# DNA 甲基化可视化 -- 工具目录

## 概述

DNA 甲基化数据可视化涵盖全局甲基化分布、染色体水平甲基化谱、DMR 热图、单基因甲基化模式、以及甲基化与表达的整合展示。

## 推荐工具

### 1. 全局甲基化水平分布

```r
library(ggplot2)
library(reshape2)

# 准备数据
meth_dist <- data.frame(
  Sample = rep(c("Ctrl_1", "Ctrl_2", "Trt_1", "Trt_2"),
               each = nrow(meth_matrix)),
  CG = c(as.matrix(meth_matrix[, 1:4])),
  Context = "CG"
)

# 密度图
ggplot(meth_dist, aes(x = CG, fill = Sample)) +
  geom_density(alpha = 0.5) +
  labs(title = "CG Methylation Level Distribution",
       x = "Methylation Level (%)",
       y = "Density") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

# 跨上下文的分布
ggplot(all_contexts, aes(x = Methylation_Level,
  fill = Context)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Sample, ncol = 2, scales = "free_y") +
  labs(title = "Methylation Levels by Context",
       x = "Methylation Level (%)") +
  theme_minimal()
```

### 2. 染色体甲基化谱

```bash
# 生成 100kb 窗口的平均甲基化
bedtools makewindows -g chrom.sizes -w 100000 \
  > genome_100kb.bed

bedtools map -a genome_100kb.bed \
  -b sample_CG_methylation.bedGraph \
  -c 4 -o mean > sample_CG_100kb.bedGraph
```

```r
# 染色体范围的甲基化谱
library(circlize)

chrom_meth <- read.table("sample_CG_100kb.bedGraph",
  col.names = c("chr", "start", "end", "methylation"))

# Circos plot
circos.initializeWithIdeogram(species = "ath")

circos.genomicTrack(chrom_meth,
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value,
      col = "steelblue", lwd = 0.5)
  },
  track.height = 0.1)
```

### 3. DMR 可视化

#### 单 DMR 甲基化信号

```r
library(DSS)

# 展示单个 DMR 附近的甲基化模式
showOneDMR(dmrs[1, ], BSobj,
  ext = 1000,   # 两侧延伸 1kb
  main = paste("DMR:", dmrs[1, "chr"],
    dmrs[1, "start"], "-", dmrs[1, "end"]))
```

#### DMR 热图

```r
library(pheatmap)

pheatmap(dmr_matrix,
  color = colorRampPalette(
    c("blue", "white", "red"))(100),
  annotation_col = sample_annotation,
  show_rownames = FALSE,
  cluster_cols = TRUE,
  main = "DMR Methylation Heatmap",
  fontsize = 8)
```

#### DMR 基因组注释分布

```r
# 注释 DMR 的基因组位置
library(ChIPseeker)

dmr_anno <- annotatePeak(makeGRangesFromDataFrame(dmrs),
  TxDb = txdb,
  tssRegion = c(-2000, 500))

# 饼图
plotAnnoPie(dmr_anno, main = "DMR Genomic Distribution")

# 到 TSS 的距离
plotDistToTSS(dmr_anno,
  title = "DMR Distribution Relative to TSS")
```

### 4. 甲基化上下文三元图

```r
# 三元图展示各上下文的甲基化比例
library(ggtern)

ternary_data <- data.frame(
  CG = cg_levels,
  CHG = chg_levels,
  CHH = chh_levels,
  Condition = condition
)

ggtern(ternary_data,
  aes(x = CG, y = CHG, z = CHH, color = Condition)) +
  geom_point(size = 2, alpha = 0.6) +
  theme_showgrid() +
  labs(title = "Methylation Context Composition",
       x = "CG%", y = "CHG%", z = "CHH%")
```

### 5. 甲基化在基因上的 Metagene Plot

```bash
# deepTools 计算基因区域的甲基化信号
computeMatrix scale-regions \
  -S sample_CG_methylation.bw \
  -R genes.bed \
  -b 2000 -a 2000 \
  --regionBodyLength 3000 \
  --binSize 50 \
  -o matrix_genes.gz

plotProfile -m matrix_genes.gz \
  -o profile_genes_methylation.pdf \
  --perGroup \
  --plotTitle "CG Methylation around Gene Body" \
  --yAxisLabel "Methylation Level (%)"

plotHeatmap -m matrix_genes.gz \
  -o heatmap_genes_methylation.pdf \
  --colorMap RdBu \
  --whatToShow "plot, heatmap and colorbar"
```

### 6. 不同上下文的比较

```r
# 箱线图展示不同上下文甲基化差异
ggplot(context_compare, aes(x = Context,
  y = Methylation_Level, fill = Condition)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(title = "Methylation Level by Context and Condition",
       x = "Sequence Context",
       y = "Methylation Level (%)") +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "salmon"))
```

### 7. 甲基化-表达关联图

```r
# 散点图 + 拟合线
ggplot(meth_expr_data, aes(x = Methylation,
  y = Expression_TPM)) +
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(method = "loess", color = "red") +
  facet_wrap(~Context) +
  labs(title = "Methylation vs Expression",
       x = "Promoter Methylation Level (%)",
       y = "Gene Expression (TPM)") +
  theme_minimal()
```

### 8. 基因组浏览器 Track

```bash
# 生成 bigWig tracks 用于 IGV/JBrowse
# CG track
bedGraphToBigWig sample_CG.bedGraph \
  chrom.sizes sample_CG.bw

# CHG track
bedGraphToBigWig sample_CHG.bedGraph \
  chrom.sizes sample_CHG.bw

# CHH track
bedGraphToBigWig sample_CHH.bedGraph \
  chrom.sizes sample_CHH.bw
```

---

## 植物特有可视化考虑

1. **三上下文分离**: 植物甲基化必须在 CG、CHG、CHH 三者之间分别展示
2. **染色体命名**: 植物染色体名如 chr1 或 Chr01 需要在可视化前统一
3. **TE 注释叠加**: 在 methylome 浏览器中叠加 TE 注释帮助解读 CHH 甲基化高峰
4. **着丝粒区域**: 着丝粒的高甲基化区域通常形成宽峰
5. **gbM 基因体甲基化**: 使用 metagene 图展示 gbM+/gbM- 基因的甲基化差异
6. **大基因组降采样**: 玉米/小麦的 1bp 分辨率甲基化可视化使用 200-500bp binning
