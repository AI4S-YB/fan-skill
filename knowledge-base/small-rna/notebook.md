# 植物小 RNA 测序分析 — 分析笔记本

## 分析概览

本笔记本提供植物小 RNA 测序数据的完整分析流程，从原始数据质控到 miRNA 功能注释。

---

## 1. 数据预处理

### 1.1 原始数据质量评估

```bash
# 使用 FastQC 评估原始测序数据质量
fastqc raw/*.fastq.gz -o qc/raw/
multiqc qc/raw/ -o qc/raw_report/
```

**检查要点**:
- Per base quality score: 所有碱基位置 Q > 20
- Adapter content: 确认接头污染程度
- Sequence length distribution: 小 RNA 文库预期为 18-30nt

### 1.2 接头去除

```bash
# Trim Galore 或 Cutadapt 去除 3' 接头
cutadapt -a TGGAATTCTCGGGTGCCAAGG \
  -m 18 -M 30 \
  --discard-untrimmed \
  -o trimmed/sample1_trimmed.fq.gz \
  raw/sample1.fq.gz
```

**植物小 RNA 接头**:
- Illumina TruSeq: TGGAATTCTCGGGTGCCAAGG
- 其他: 根据建库试剂盒选择

### 1.3 长度筛选与质量控制

```bash
# 去除低质量 reads (Q < 20)
fastq_quality_filter -q 20 -p 80 -i trimmed/sample1_trimmed.fq.gz \
  -o filtered/sample1_filtered.fq.gz

# 转换为 FASTA 格式
seqkit fq2fa filtered/sample1_filtered.fq.gz \
  -o filtered/sample1_filtered.fa
```

### 1.4 去除非 miRNA 序列

```bash
# 比对到 rRNA/tRNA/snRNA 数据库，去除污染
bowtie -v 0 -k 1 --un filtered/sample1_noncoding.fa \
  rfam_rrna_trna_index filtered/sample1_filtered.fa \
  > /dev/null
```

---

## 2. miRNA 鉴定

### 2.1 已知 miRNA 鉴定 (模式物种)

```bash
# miRDeep2 使用 miRBase 参考
mapper.pl config.txt -d -e -h -m \
  -p genome_index \
  -s filtered/sample1_noncoding.fa \
  -t filtered/sample1_noncoding.arf

miRDeep2.pl filtered/sample1_noncoding.fa \
  genome.fa \
  filtered/sample1_noncoding.arf \
  mirbase_mature.fa \
  mirbase_precursor.fa \
  mirbase_others.fa \
  -t Osa 2> report.log
```

### 2.2 新 miRNA 预测 (非模式物种)

```bash
# 不提供已知 miRNA 参考，仅使用基因组
miRDeep2.pl filtered/sample1_noncoding.fa \
  genome.fa \
  filtered/sample1_noncoding.arf \
  none none none \
  -t Novel 2> report.log
```

**结果筛选标准**:
- miRDeep2 score > 0 (推荐 > 4)
- 显著 Randfold p-value (p < 0.05)
- 成熟的 miRNA 序列不在 rRNA/tRNA 区域

---

## 3. miRNA 靶基因预测

### 3.1 psRNATarget

```bash
# 使用 psRNATarget 在线服务器或本地版本
# 输入: miRNA.fasta, transcriptome.fasta
```

**参数设置**:
- Expectation ≤ 5: 默认阈值
- Maximum mismatches: ≤ 4 (植物 miRNA-mRNA 高度互补)
- 翻译抑制范围: 9-11nt 位点 (剪切位点)

### 3.2 降解组分析 (靶基因验证)

```bash
# CleaveLand4 鉴定剪切位点
cleaveland4.pl \
  -e degradome.fa \
  -u "T" \
  -t transcriptome.fa \
  -p miRNA.fa \
  -o cleaveland_output/
```

**剪切位点分类**:
- Category 0: 只有一个峰，>1 RPM
- Category 1: 最高峰且 >1 RPM
- Category 2: 不是最高峰，但 >1 RPM
- Category 3: >1 RPM 但不在第 10 位
- Category 4: <1 RPM

---

## 4. 差异表达分析

### 4.1 表达定量

```bash
# 使用 featureCounts 或 HTSeq 定量 miRNA 表达
featureCounts -a mirna.gff3 -o counts.txt \
  -t miRNA -g ID \
  alignments/*.bam
```

### 4.2 DESeq2 差异分析

```r
library(DESeq2)

# 读取 count 矩阵
counts <- read.table("counts.txt", header = TRUE, row.names = 1)

# 实验设计
colData <- data.frame(
  row.names = colnames(counts),
  condition = c("control", "control", "control", "treatment", "treatment", "treatment")
)

# DESeq2 分析
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData,
                               design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treatment", "control"))
res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# 导出结果
write.csv(as.data.frame(res_sig), "diff_mirna_results.csv")
```

---

## 5. 可视化

### 5.1 差异 miRNA 火山图

```r
library(ggplot2)
library(ggrepel)

res$sig <- "NS"
res$sig[res$padj < 0.05 & res$log2FoldChange > 1] <- "Up"
res$sig[res$padj < 0.05 & res$log2FoldChange < -1] <- "Down"

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down" = "blue", "NS" = "grey", "Up" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "差异 miRNA 火山图", x = "log2 Fold Change", y = "-log10(padj)") +
  theme_bw()
```

### 5.2 差异 miRNA 热图

```r
library(pheatmap)

# 选取差异 miRNA 并绘制热图
top_mirnas <- rownames(res_sig)[1:50]
pheatmap(norm_counts[top_mirnas, ],
         scale = "row",
         annotation_col = colData,
         main = "差异 miRNA 表达热图")
```

---

## 常见问题

### Q: 比对率低怎么办?
A: 检查是否去除了接头，确认参考基因组版本是否正确。植物小 RNA 数据比对率通常为 60-80%。

### Q: 预测的 miRNA 中包含大量 rRNA 片段怎么办?
A: 确保在预处理阶段使用了完整的 rRNA/tRNA 数据库进行过滤。

### Q: 降解组分析无显著结果怎么办?
A: 考虑放宽 RNA 丰度阈值 (RPM cutoff)，或检查 miRNA-mRNA 互补性。

---

## 植物小 RNA 的特有生物学

### 植物 miRNA 与动物 miRNA 的关键差异

植物 miRNA 分析不能照搬动物 miRNA 流程，两者有本质区别：

| 特性 | 植物 miRNA | 动物 miRNA |
|------|-----------|-----------|
| 靶基因配对 | 近乎完全互补（<4错配） | 主要通过seed区域（6-8nt） |
| 作用机制 | 以mRNA剪切为主 | 以翻译抑制为主 |
| 基因组位置 | 独立转录单位（MIR基因） | 多位于内含子中 |
| 保守性 | 家族保守（miR156/172/159等古老家族） | 家族保守但序列变异大 |
| 长度 | 21nt为主 | 21-23nt |
| 甲基化 | 3'端2'-O-甲基化（HEN1介导） | 无（Dicer产物无修饰） |

**关键影响**：
- 靶基因预测方法完全不同：植物用 psRNATarget（高互补性模型），动物用 TargetScan/miRanda（seed配对模型）。**绝不能混用**。
- 植物 miRNA-mRNA 的高互补性意味着更少的靶基因（每个 miRNA 通常 1-10 个靶标），但每个靶标的调控效应更强。

### 植物特有的小 RNA 类别

植物中有多种动物中没有或功能不同的 sRNA 类别：

**phasiRNA（phased secondary siRNA）**：
- 产生机制：miRNA 剪切特定非编码转录本（PHAS/TAS）后，RDR6 合成 dsRNA，由 DCL4 按 21nt 相位切割
- 功能：生殖发育（花药/花粉）、逆境响应
- 鉴定工具：PhaseTank、ShortStack
- 植物特有且生物学意义重要，不应忽视

**ta-siRNA（trans-acting siRNA）**：
- 是 phasiRNA 的一个子类，其前体为 TAS 基因家族
- TAS3 靶向多个 ARF（auxin response factor）转录本，调控生长素信号转导
- 几乎所有陆生植物都保守

**24nt siRNA（heterochromatic siRNA, hc-siRNA）**：
- 植物 sRNA 文库中最丰富的组分（占总 reads 30-60%）
- 源自重复序列和转座子区域
- 通过 RdDM（RNA-directed DNA Methylation）途径引导 DNA 甲基化
- 对维持基因组稳定性和 TE 沉默至关重要
- 在常规 miRNA 分析中通常被过滤（因长度 > 22nt），但如果你关注表观遗传调控，24nt siRNA 是重要数据源

### 保守 miRNA 家族 vs 物种特异性 miRNA

植物 miRNA 分为两大类：

**保守 miRNA 家族**：
- miR156/157（靶向 SPL 转录因子，调控开花时间）
- miR159/319（靶向 MYB/TCP，调控叶片发育和激素应答）
- miR160（靶向 ARF，调控生长素信号）
- miR164（靶向 NAC，调控器官边界和侧根发育）
- miR165/166（靶向 HD-ZIP III，调控维管束和叶极性发育）
- miR172（靶向 AP2，调控开花和花器官发育）
- miR396（靶向 GRF，调控细胞增殖）
- 这些 miRNA 在几乎所有陆生植物中保守 —— 如果这些家族在你的非模式物种中"不见了"，多半是注释问题而非真的缺失

**物种 / 谱系特异性 miRNA**：
- 通常为近期演化产生的 miRNA
- 可能没有 miRBase 条目（非模式物种）
- 靶基因富集在某些特定生物学过程（如次生代谢调控、共生信号）
- 表达水平通常低于保守 miRNA

### 多倍体植物 miRNA 分析

多倍体中同一 miRNA 家族可能有来自不同亚基因组的多个 copies：

- **区分 homeolog copies**：如果已有亚基因组组装，将 sRNA reads 分别比对到各亚基因组
- **表达偏倚**（homeolog expression bias）：一个亚基因组的 miRNA copy 表达量可能显著高于另一个。这在多倍体演化中有重要生物学意义
- **靶基因同步分化**：不同亚基因组的 miRNA copies 可能靶向不同的 homeolog 靶基因（subfunctionalization）
- **保守性判断**：在多倍体中，保守 miRNA 的 copies 数大致等于 ploidy 水平（如四倍体棉花中每个保守 miRNA 家族约 4 个 copies）

### 作物特异性 miRNA 资源

| 资源 | 适用作物 | 内容 |
|------|---------|------|
| PMRD（Plant miRNA Database） | 所有植物 | 覆盖 120+ 植物物种的 miRNA 序列和靶基因 |
| PmiREN（Plant miRNA ENcyclopedia） | 所有植物 | 6000+ 植物物种的基因组范围 miRNA 注释 |
| miRBase | 模式植物 | 高置信度 curated miRNA（拟南芥 ~400，水稻 ~600） |
| sRNAanno | 所有植物 | 基于基因组比对的注释，适合非模式物种 |
| mirtronPred | 所有植物 | 预测 mirtron 类型 miRNA（位于内含子中） |

对于非模式作物（如茶 Camellia sinensis、木薯 Manihot esculenta、小米 Setaria italica），PmiREN 比 miRBase 提供更完整的 miRNA 注释参考。在发表时，使用 PmiREN 或 sRNAanno 的 miRNA 家族编号（而非自己随意命名）以确保交叉引用。

---

## 参考

- miRDeep2: https://github.com/rajewsky-lab/mirdeep2
- psRNATarget: http://plantgrn.noble.org/psRNATarget/
- CleaveLand: https://github.com/MikeAxtell/CleaveLand4
- DESeq2: https://bioconductor.org/packages/DESeq2/
