# DNA 甲基化数据分析 -- 分析笔记本

## 分析概览

本笔记本提供植物 DNA 甲基化数据的完整分析流程，涵盖全基因组亚硫酸盐测序（WGBS）和简化代表性亚硫酸盐测序（RRBS）的数据处理、甲基化水平定量、差异甲基化区域（DMR）检测及功能注释。

---

## 1. 数据预处理

### 1.1 原始数据质控

```bash
# FastQC 评估原始数据质量
fastqc raw/*.fastq.gz -o qc/raw/
multiqc qc/raw/ -o qc/raw_report/
```

**检查要点**:
- Per base quality score: Q > 20
- GC content: WGBS/RRBS 的 GC 分布不同于常规 DNA-seq（亚硫酸盐转化影响核酸组成）
- Bismark 后续分析可耐受适度质量下降

### 1.2 接头与质量修剪

```bash
# Trim Galore 专门用于 bisulfite 数据的修剪
trim_galore --paired --clip_R1 10 --clip_R2 10 \
  --three_prime_clip_R1 10 --three_prime_clip_R2 10 \
  --fastqc \
  sample_R1.fq.gz sample_R2.fq.gz \
  -o trimmed/
```

**注意事项**:
- WGBS 数据建议去接头而非质量修剪（亚硫酸盐转化引入的 C->T 转换可能被误判为低质量）
- 默认 Trim Galore 质量阈值 Q20
- RRBS 数据需要 `--rrbs` 参数处理酶切末端

---

## 2. 比对与甲基化提取

### 2.1 Bismark 比对 (WGBS)

```bash
# 基因组 bisulfite 转换索引
bismark_genome_preparation --path_to_aligner bowtie2 \
  --bowtie2 genome_dir/

# Bismark 比对
bismark --bowtie2 -p 4 --non_directional \
  --genome genome_dir/ \
  -1 sample_R1_trimmed.fq.gz \
  -2 sample_R2_trimmed.fq.gz \
  --output_dir bismark_output/

# 去重复
deduplicate_bismark -p --bam \
  bismark_output/sample_bismark_bt2_pe.bam \
  --output_dir bismark_output/
```

### 2.2 甲基化提取

```bash
# 提取全基因组甲基化信息
bismark_methylation_extractor --paired-end --no_overlap \
  --bedGraph --CX_context \
  --cytosine_report --genome_folder genome_dir/ \
  bismark_output/sample_deduplicated.bam \
  --output_dir methylation_output/

# 生成 cytosine report
coverage2cytosine --genome_folder genome_dir/ \
  -o sample_CX_report.txt \
  methylation_output/sample_deduplicated.bismark.cov.gz
```

### 2.3 Bismark 比对 (RRBS)

```bash
# RRBS 模式比对 (添加 --rrbs 参数)
bismark --bowtie2 -p 4 \
  --genome genome_dir/ \
  -1 sample_R1_trimmed.fq.gz \
  -2 sample_R2_trimmed.fq.gz \
  --rrbs \
  --output_dir bismark_output/
```

### 2.4 植物特有甲基化上下文

植物 DNA 甲基化发生在三种序列上下文中：

- **CG (CpG)**: 对称甲基化，由 MET1 维持
- **CHG**: 对称甲基化，由 CMT3 维持
- **CHH**: 非对称甲基化，由 DRM2 维持（RdDM 途径）

```bash
# 分离不同上下文的甲基化水平
# CG 上下文
grep "CG_" sample_CX_report.txt > sample_CG_report.txt

# CHG 上下文
grep "CHG_" sample_CX_report.txt > sample_CHG_report.txt

# CHH 上下文
grep "CHH_" sample_CX_report.txt > sample_CHH_report.txt
```

**各上下文典型甲基化水平**:
| 上下文 | 拟南芥叶片 | 水稻叶片 | 玉米叶片 |
|--------|-----------|---------|---------|
| CG | 24-30% | 30-45% | 50-60% |
| CHG | 6-10% | 10-15% | 20-30% |
| CHH | 2-5% | 3-8% | 2-5% |

---

## 3. 甲基化水平定量

### 3.1 单碱基分辨率甲基化水平

```r
library(methylKit)

# 读取 cytosine report
file.list <- list("sample1_CG.txt", "sample2_CG.txt",
                  "sample3_CG.txt", "sample4_CG.txt")

myobj <- methRead(file.list,
                  sample.id = list("sample1", "sample2", "sample3", "sample4"),
                  assembly = "genome",
                  treatment = c(0, 1, 0, 1),
                  context = "CpG")
```

### 3.2 甲基化水平统计

```r
# 全局甲基化水平
getMethylationStats(myobj[[1]], plot = TRUE, both.strands = FALSE)

# 覆盖率统计
getCoverageStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
```

---

## 4. DMR 检测

### 4.1 DSS DMR Calling

```r
library(DSS)

# 读取各样本的甲基化文件
sample1 <- read.table("sample1_CG_report.txt", header = TRUE)
sample2 <- read.table("sample2_CG_report.txt", header = TRUE)
sample3 <- read.table("sample3_CG_report.txt", header = TRUE)
sample4 <- read.table("sample4_CG_report.txt", header = TRUE)

# 构建 BSseq 对象
BSobj <- makeBSseqData(
  list(sample1, sample2, sample3, sample4),
  c("sample1", "sample2", "sample3", "sample4")
)

# DML 检验
dmlTest <- DMLtest(BSobj,
                   group1 = c("sample1", "sample3"),
                   group2 = c("sample2", "sample4"),
                   smoothing = TRUE)

# 提取 DML
dmls <- callDML(dmlTest, p.threshold = 1e-5)

# 检测 DMR
dmrs <- callDMR(dmlTest, p.threshold = 1e-5,
                minlen = 50, minCG = 3,
                dis.merge = 100)

# 查看结果
head(dmrs)
write.csv(dmrs, "dmr_results.csv")
```

**DSS 参数说明**:
- `p.threshold`: DML/DMR 的显著性阈值
- `minlen`: DMR 最小长度（bp）
- `minCG`: DMR 最少包含的 CpG 位点数
- `dis.merge`: 相邻 DML 合并距离

---

## 5. 甲基化-表达关联分析

### 5.1 基因启动子甲基化与表达关联

```r
# 读取表达数据
expression <- read.csv("gene_expression.csv", row.names = 1)

# 注释甲基化位点到基因启动子区域
library(annotatr)
annots <- build_annotations(genome = "tair10",
                            annotations = "tair10_basicgenes")

# 关联分析
correlation_results <- data.frame()
for (gene in intersect(rownames(expression), gene_promoter_meth$gene)) {
  meth_levels <- gene_promoter_meth[gene_promoter_meth$gene == gene, "methylation"]
  expr_levels <- expression[gene, ]
  cor_test <- cor.test(meth_levels, expr_levels, method = "spearman")
  # 存储结果
}
```

### 5.2 转录因子结合与甲基化

```r
# 低甲基化区域 (Hypo-DMR) 的 TF 结合 motif 富集
# 使用 HOMER 或 MEME 分析 DMR 区域的 motif 富集
```

**植物甲基化调控的特殊性**:
- DNA 甲基化在启动子区域通常与转录抑制相关
- 基因体甲基化（gene body methylation）在植物中与中等表达水平相关
- CHH 甲基化主要由 RdDM 途径维护，与 TE 沉默相关
- CG 甲基化在基因体中的功能与动物不同（植物中 gbM 与表达稳定性相关）

---

## 6. 可视化

### 6.1 甲基化水平分布

```r
# 全局甲基化水平分布
library(ggplot2)
ggplot(methylation_data, aes(x = methylation_level, fill = context)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~sample) +
  theme_minimal() +
  labs(title = "Methylation Level Distribution",
       x = "Methylation Level (%)", y = "Density")
```

### 6.2 基因组浏览器 Track

```bash
# 生成 bigWig track
bedGraphToBigWig sample_CG.bedGraph chrom.sizes sample_CG.bw
```

### 6.3 DMR 热图

```r
library(pheatmap)
pheatmap(dmr_matrix,
         annotation_col = sample_annotation,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE,
         main = "DMR Methylation Levels")
```

---

## 常见问题

### Q: Bismark 比对率低怎么办？
A: WGBS 数据的唯一比对率通常低于常规 DNA-seq（亚硫酸盐转化降低了序列复杂度）。40-60% 是正常范围。确保使用 `--non_directional` 参数（植物全基因组甲基化是非方向性的）。

### Q: CHH 甲基化水平过低？
A: CHH 甲基化通常较低（1-5%），需要较深的测序深度（>30x）才能准确检测。RRBS 方法可能无法捕获足够的 CHH 位点。

### Q: DMR 数目过多或过少？
A: 调整 DSS 的 `p.threshold` 参数。过少尝试 1e-3，过多尝试 1e-8。同时检查生物学重复的一致性。

### Q: 如何区分 CHG 和 CHH ？
A: Bismark 的 CX_report 输出文件自动区分三种上下文（CG/CHG/CHH）。CHG 中的 H 代表 A、T 或 C（非 G），CHH 中的两个 H 都不能是 G。

---

## 参考

- Bismark: https://github.com/FelixKrueger/Bismark
- DSS: https://bioconductor.org/packages/DSS/
- methylKit: https://bioconductor.org/packages/methylKit/
- Trim Galore: https://github.com/FelixKrueger/TrimGalore
- Plant Methylation DB: http://epigenome.genetics.uga.edu/PlantMethylDB/
