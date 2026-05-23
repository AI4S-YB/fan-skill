# 基因组组装可视化 -- 工具目录

## 概述

基因组组装可视化涵盖组装图结构、contig 连续性、Hi-C contact map、基因组完整性概览和 GC 含量/覆盖度散点图。

## 推荐工具

### 1. Bandage -- 组装图可视化

**描述**: Bandage 是组装图（GFA 格式）的交互式可视化工具，可显示 contig 连接关系、覆盖度、和重复结构。

**安装**:
```bash
conda install -c bioconda bandage
```

**命令行用法**:
```bash
# 生成高质量 PNG
Bandage image assembly_graph.gfa assembly_graph.png \
  --height 3000 --width 4000 \
  --lengths --colour random

# 交互式界面
Bandage load assembly_graph.gfa
```

**在 Bandage 中**:
- 每条边 = 一个 contig (连续序列)
- 节点间的线 = overlapping 或 bridging 关系
- 颜色和宽度反映覆盖度和长度

### 2. Hi-C Contact Map (HiCExplorer)

```bash
# 构建 Hi-C 交互矩阵
hicBuildMatrix -s hic_sorted.bam \
  --binSize 50000 \
  --inputBufferSize 4000 \
  -o hic_matrix.h5

# 绘制 contact map
hicPlotMatrix -m hic_matrix.h5 \
  -o hic_contact_map.png \
  --log1p \
  --dpi 300 \
  --title "Assembly Hi-C Contact Map" \
  --colorMap RdYlBu_r
```

**解读 Hi-C Contact Map**:
- **对角线强信号**: 染色体内部连续性良好
- **清晰的正方形块**: 染色体间的边界分明
- **"plume" 或 "haze"**: 可能的 misassembly 或真实染色体互作

### 3. Coverage-GC 散点图

```bash
# 使用 blobtools 生成 Coverage-GC 图
# 可用于检测污染 (叶绿体、线粒体、细菌)
blobtools create -i assembly.fasta \
  -b coverage.bam \
  -t coverage.tsv \
  -o blobtools_output/

blobtools view -i blobtools_output.blobDB.json \
  -o blobtools_output/

blobtools plot -i blobtools_output.blobDB.json \
  -o blobtools_output/
```

### 4. 累积长度图 (Nx 曲线)

```r
library(ggplot2)

# 读取 contig 长度
contigs <- read.table("contig_lengths.txt",
  col.names = "length")
contigs <- contigs[order(contigs$length, decreasing = TRUE), ]

# 计算 Nx
contigs$cum_length <- cumsum(contigs$length)
contigs$cum_frac <- contigs$cum_length / sum(contigs$length)

# Nx 图
ggplot(contigs, aes(x = 1:nrow(contigs), y = cum_frac)) +
  geom_step(size = 1, color = "steelblue") +
  geom_hline(yintercept = 0.5, linetype = "dashed",
    color = "red", alpha = 0.5) +
  labs(title = "Assembly Nx Curve",
       x = "Contig Rank (#)",
       y = "Cumulative Fraction of Assembly") +
  theme_minimal()
```

### 5. BUSCO 总结图

```bash
# BUSCO 自动生成可视化
# output: busco_output/run_*/busco_figure.png

# 或使用 R 自定义
```

```r
# BUSCO 条形图
library(ggplot2)

busco_results <- data.frame(
  Category = c("Complete (S)", "Complete (D)",
    "Fragmented", "Missing"),
  Percentage = c(82.1, 8.1, 4.3, 5.5)
)

ggplot(busco_results, aes(x = Category, y = Percentage,
  fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("green4", "green3",
    "orange", "red3")) +
  labs(title = "BUSCO Assessment Results",
       y = "Percentage (%)") +
  theme_minimal() +
  theme(legend.position = "none")
```

### 6. Merqury Spectra Plot

```bash
# Merqury 自动生成的频谱图
# 文件: output_prefix.spectra-cn.png

# 频谱图解读:
# - 峰 1 (1x): 杂合位点/测序错误
# - 峰 2 (2x): 二倍体基因组
# - 峰 0.5 (0.5x): 低覆盖重复/错误
```

### 7. 各染色体可视化

```r
library(karyoploteR)

# 染色体核型图
karyoplot <- plotKaryotype(genome = "tair10",
  chromosomes = c("1", "2", "3", "4", "5"))

# 添加 BUSCO 位置
kpPlotMarkers(karyoplot,
  chr = busco_positions$chr,
  x = busco_positions$pos,
  labels = busco_positions$busco_id,
  marker.parts = c(1, 0, 0),
  cex = 0.3)

# 添加 gap (N) 位置
kpPlotRegions(karyoplot,
  data = gap_regions,
  col = "red")
```

### 8. 基因组组装仪表板

```bash
# 使用 QUAST 生成完整的组装报告
quast assembly.fasta \
  -r reference_genome.fasta \
  -g genes.gff \
  --eukaryote \
  --large \
  -o quast_output/

# QUAST 输出包括:
# - HTML 报告
# - 累积长度图
# - missassembly 分析
# - 基因区域覆盖率
```

---

## 植物组装可视化特殊考虑

1. **大型基因组的降采样**: 玉米 (2.4Gb) 或小麦 (17Gb) 的 Hi-C contact map 需用 100-500kb bin 以避免噪音
2. **亚基因组的颜色编码**: 异源多倍体（如面包小麦 6x）的染色体应按亚基因组分组着色（A 基因组 = 红, B = 蓝, D = 绿）
3. **细胞器基因组**: 在 Coverage-GC 图中，叶绿体和线粒体的高覆盖度点通常偏离主簇，可用于鉴定细胞器 contig
4. **植物特有的图结构**: 着丝粒和 NOR (Nucleolar Organizer Region) 在组装图中可能形成"菊花"或"星形"结构
5. **多重度标记**: 在 Bandage 中可通过节点颜色标注重复度，帮助识别 collapsed repeats

---

## 输出建议

**论文用图**:
- 分辨率: 300-600 DPI
- 格式: PDF (矢量) 或 TIFF (光栅)
- 配色: 无障碍色板 (colorbrewer)

**交互式探索**:
- 格式: HTML (QUAST), 交互式 IGV, JBrowse2
