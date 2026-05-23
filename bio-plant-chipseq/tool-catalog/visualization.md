# ChIP-seq 可视化 — 工具目录

## 概述

ChIP-seq 数据的可视化对于结果展示和数据质量评估至关重要。本目录涵盖常用的可视化类型和推荐工具。

## 可视化类型与工具

### 1. 基因组浏览器可视化

**工具**: IGV (Integrative Genomics Viewer)

**文件格式**:
- **bigWig**: 连续信号 track (覆盖度)
- **BED**: 离散 peak track

**生成 bigWig 文件 (deepTools)**:
```bash
bamCoverage -b sample_sorted.bam -o sample.bw \
  --normalizeUsing RPKM \
  --binSize 10
```

### 2. Peak 热图与 Meta-gene Profile

**工具**: deepTools (plotHeatmap, plotProfile)

**描述**: 展示多个样品在 TSS、peak 中心或指定区域的信号分布。

**生成 TSS 附近信号矩阵**:
```bash
computeMatrix reference-point -S sample1.bw sample2.bw \
  -R genes.bed \
  --referencePoint TSS \
  -a 3000 -b 3000 \
  -o matrix_TSS.gz
```

**绘制热图**:
```bash
plotHeatmap -m matrix_TSS.gz \
  -o heatmap_TSS.pdf \
  --colorMap Blues
```

**绘制 Profile**:
```bash
plotProfile -m matrix_TSS.gz \
  -o profile_TSS.pdf \
  --perGroup
```

### 3. Peak 分布饼图

**工具**: ChIPseeker (R)

**描述**: 展示 peaks 在基因组功能区域的分布比例。

**R 代码**:
```r
library(ChIPseeker)
plotAnnoPie(peak_anno)
plotAnnoBar(peak_anno)
```

### 4. Fingerprint Plot (指纹图)

**工具**: deepTools (plotFingerprint)

**描述**: 评估 ChIP-seq 样本的信号富集度，展示 reads 在基因组上的累积分布。

```bash
plotFingerprint -b sample_chip.bam sample_input.bam \
  --labels ChIP Input \
  -o fingerprint.pdf
```

### 5. 差异结合可视化

#### Volcano Plot
```r
library(ggplot2)
ggplot(db_results, aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() +
  labs(title = "差异结合位点火山图")
```

#### MA Plot
**工具**: DiffBind (dba.plotMA)

#### PCA Plot
**工具**: DiffBind (dba.plotPCA) 或 deepTools (plotPCA)

### 6. Motif Logo 可视化

**工具**: MEME Suite

**描述**: 展示富集的 DNA motif 序列 logos。

---

## 植物特殊注意事项

1. **插入/缺失区域**: 植物基因组与参考基因组的差异可能导致 reads 堆积异常
2. **叶绿体/线粒体基因组**: 排除来自细胞器基因组的 reads
3. **大染色体可视化**: 对于大基因组物种（玉米、小麦），可能需要按染色体分别可视化
