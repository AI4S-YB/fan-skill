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

## 参考

- miRDeep2: https://github.com/rajewsky-lab/mirdeep2
- psRNATarget: http://plantgrn.noble.org/psRNATarget/
- CleaveLand: https://github.com/MikeAxtell/CleaveLand4
- DESeq2: https://bioconductor.org/packages/DESeq2/
