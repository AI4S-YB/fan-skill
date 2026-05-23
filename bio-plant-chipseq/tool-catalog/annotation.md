# Peak 注释 — 工具目录

## 概述

Peak 注释将 ChIP-seq peaks 关联到附近的基因，是理解 ChIP-seq 结果生物学功能的关键步骤。

## 推荐工具

### 1. ChIPseeker (推荐)

**描述**: R/Bioconductor 包，用于 ChIP-seq peak 注释和可视化。

**关键功能**:
- Peak 到基因的注释 (promoter, 5' UTR, 3' UTR, exon, intron, intergenic)
- TSS 附近 peak 分布可视化
- 功能富集分析 (GO, KEGG)

**R 代码示例**:
```r
library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)  # 拟南芥
# 或使用 GFF/GTF 文件
# txdb <- makeTxDbFromGFF("genome.gff3")

# 读取 peaks
peaks <- readPeakFile("tf_peaks.narrowPeak")

# 注释 peaks
peak_anno <- annotatePeak(peaks,
                          TxDb = txdb,
                          tssRegion = c(-3000, 3000))

# 查看注释分布
plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno)

# 提取注释结果
anno_df <- as.data.frame(peak_anno)
```

**植物 TxDb 资源**:
| 物种 | R 包/文件 |
|------|----------|
| 拟南芥 | `TxDb.Athaliana.BioMart.plantsmart28` |
| 水稻 | `TxDb.Osativa.MSU.msu7` 或 GFF 文件 |
| 玉米 | 自定义 makeTxDbFromGFF |
| 大豆 | 自定义 makeTxDbFromGFF |
| 其他 | makeTxDbFromGFF("species.gff3") |

**注释分类**:
- **Promoter (<= 1kb)**: 启动子近端区域
- **Promoter (1-2kb)**: 启动子远端区域
- **Promoter (2-3kb)**: 启动子最远端区域
- **5' UTR**: 5‘ 非翻译区
- **3' UTR**: 3' 非翻译区
- **Exon**: 外显子
- **Intron**: 内含子
- **Downstream (<= 3kb)**: 基因下游
- **Distal Intergenic**: 基因间远距离

### 2. bedtools closest

**描述**: 命令行工具，用于快速找出距离每个 peak 最近的基因。

```bash
bedtools closest -a peaks.bed -b genes.bed -d > peak_gene_distances.txt
```

### 3. HOMER (annotatePeaks.pl)

**描述**: HOMER 的 peak 注释功能，直接输出注释表格。

**示例**:
```bash
annotatePeaks.pl peaks.bed tair10 > peak_annotation.txt
```

---

## 植物特殊注意事项

1. **基因间区 peaks**: 植物基因组中远距离调控元件较常见，注意基因间区 peaks 的功能
2. **非编码 RNA**: 植物中存在大量 lncRNA，可能需要额外注释
3. **转座子区域**: 植物基因组中转座子较多，组蛋白修饰 peaks 可能富集在 TE 区域
4. **启动子定义**: 植物 TSS 注释不如动物完善，建议使用较大的启动子窗口 (>2kb)
5. **基因密度**: 植物基因组基因密度差异大，影响 peak 注释的准确性
