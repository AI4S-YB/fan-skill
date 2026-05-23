# 甲基化水平定量 -- 工具目录

## 概述

甲基化水平定量是从 Bismark 输出的 cytosine report 中计算每个胞嘧啶位点或每个区域的甲基化水平，并进一步统计基因组特征（如基因、启动子、TE）的甲基化谱。

## 推荐工具

### 1. methylKit -- R 工具包

**描述**: methylKit 是甲基化数据分析和统计的综合性 R 包，支持多种输入格式和丰富的分析功能。

**安装**:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("methylKit")
```

### 单碱基分辨率甲基化水平

```r
library(methylKit)

# 读取各样本 cytosine report
file.list <- list(
  "sample1_CG.txt", "sample2_CG.txt",
  "sample3_CG.txt", "sample4_CG.txt"
)

myobj <- methRead(file.list,
  sample.id = list("ctrl_rep1", "ctrl_rep2",
                   "trt_rep1", "trt_rep2"),
  assembly = "genome",
  treatment = c(0, 0, 1, 1),
  context = "CpG",
  mincov = 10)   # 最少覆盖度
```

**参数说明**:
| 参数 | 说明 |
|------|------|
| `mincov` | 最少 reads 覆盖度（推荐 10x） |
| `pipeline` | 数据来源：bismark (默认) |
| `context` | CpG, CHG, CHH (分别读取) |

### 基础统计

```r
# 全局甲基化水平
for (i in 1:length(myobj)) {
  cat("Sample", myobj[[i]]@sample.id, "\n")
  getMethylationStats(myobj[[i]],
    plot = TRUE, both.strands = FALSE)
}

# 覆盖率统计
for (i in 1:length(myobj)) {
  getCoverageStats(myobj[[i]],
    plot = TRUE, both.strands = FALSE)
}
```

### 样本间甲基化水平比较

```r
# 合并所有样本的共享位点
meth <- unite(myobj, destrand = FALSE)
head(meth)

# 查看合并后的位点数
cat("Shared sites after merging:", nrow(meth), "\n")

# 各样本甲基化水平分布
par(mfrow = c(2, 2))
for (i in 1:ncol(meth)) {
  hist(meth[, i],
    breaks = 100, main = paste("Sample", i),
    xlab = "Methylation level (%)",
    col = "steelblue")
}
```

### 按基因组特征定量

```r
# 读取基因组注释（GFF/BED）
library(genomation)

# 基因启动子区域 (TSS +/- 1kb)
promoters <- readGeneric("promoters_1kb.bed",
  header = FALSE,
  keep.all.metadata = TRUE)

# 计算启动子区域的甲基化水平
meth_promoter <- regionCounts(myobj, promoters)

# 基因体甲基化 (gene body methylation, gbM)
gene_bodies <- readGeneric("genes.bed")
meth_genebody <- regionCounts(myobj, gene_bodies)

# 转座子 (TE) 甲基化
tes <- readGeneric("transposons.bed")
meth_te <- regionCounts(myobj, tes)
```

### 2. 自定义脚本 -- 甲基化水平矩阵

使用 Bismark 的 cytosine report 生成全基因组甲基化矩阵：

```bash
# 提取甲基化水平，生成 BED-like 格式
awk 'BEGIN{OFS="\t"}
  NR>1 {
    methylated = $4;
    unmethylated = $5;
    total = methylated + unmethylated;
    if(total >= 5) {
      level = methylated / total * 100;
      print $1, $2, $3, level, methylated, unmethylated;
    }
  }' sample_CG_report.txt > sample_CG_methylation.bedGraph
```

### 3. 滑动窗口甲基化水平

```bash
# 使用 bedtools 计算非重叠窗口的平均甲基化
bedtools makewindows -g chrom.sizes -w 10000 \
  > genome_10kb_windows.bed

bedtools map -a genome_10kb_windows.bed \
  -b sample_CG_methylation.bedGraph \
  -c 4 -o mean > sample_CG_10kb_mean.bedGraph
```

---

## 各上下文的甲基化水平统计

### 全局甲基化水平 (植物参考值)

```r
# 计算各上下文全局甲基化
for (context in c("CG", "CHG", "CHH")) {
  obj <- methRead(list("sample_report.txt"),
    context = context, mincov = 5)
  mean_meth <- mean(obj@methylation)
  cat(context, "global methylation:",
    round(mean_meth, 1), "%\n")
}
```

**植物参考值（叶片组织）**:

| 物种 | CG (%) | CHG (%) | CHH (%) |
|------|--------|---------|---------|
| 拟南芥 | 24-30 | 6-10 | 2-5 |
| 水稻 | 30-45 | 10-15 | 3-8 |
| 玉米 | 50-60 | 20-30 | 2-5 |
| 大豆 | 40-55 | 15-25 | 3-6 |
| 番茄 | 50-65 | 20-30 | 5-10 |

### Gene Body Methylation (gbM)

```r
# 计算基因体的甲基化
# gbM 是植物的独特特征，某些基因的基因体被甲基化

gbM_score <- apply(gene_meth_matrix, 1, mean)
gbM_genes <- names(gbM_score[gbM_score > 20])

cat("Number of gbM genes:", length(gbM_genes), "\n")
cat("Fraction of total genes:",
  round(length(gbM_genes) / nrow(gene_meth_matrix) * 100, 1),
  "%\n")
```

**gbM 特征**:
- 主要存在于 CG 上下文中
- 与外显子位置相关（不覆盖内含子）
- 排除转录起始位点（TSS）和转录终止位点（TES）
- 功能上与中等表达水平、低表达噪声相关

---

## 植物特有注意事项

1. **低覆盖位点过滤**: CHH 甲基化水平低，需要 `mincov >= 10` 而非 5
2. **CG 位点分布**: 植物 CG 密度远低于动物，尤其在基因间区
3. **着丝粒/端粒**: 着丝粒附近的甲基化水平可能与基因组平均水平显著不同
4. **细胞器 DNA**: 线粒体和叶绿体 DNA 几乎不甲基化，可作为非甲基化对照
5. **CHH 与 24nt siRNA**: CHH 甲基化主要由 24-nt siRNA 引导（RdDM 途径），在 TE 附近和年轻插入序列中显著升高
