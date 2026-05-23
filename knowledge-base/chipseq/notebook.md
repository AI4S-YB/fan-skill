# 植物 ChIP-seq 数据分析 — 分析笔记本

## 分析概览

本笔记本提供植物 ChIP-seq 数据的完整分析流程，涵盖转录因子和组蛋白修饰 ChIP-seq 的数据处理、peak calling、差异分析和功能注释。

---

## 1. 数据预处理

### 1.1 原始数据质控

```bash
# 使用 FastQC 评估原始数据质量
fastqc raw/*.fastq.gz -o qc/raw/
multiqc qc/raw/ -o qc/raw_report/
```

**检查要点**:
- Per base quality score: Q > 20
- GC content: 与参考基因组 GC 含量接近
- Adapter content: 确认无接头污染

### 1.2 Reads 比对到参考基因组

```bash
# 使用 Bowtie2 进行比对 (植物推荐 --very-sensitive)
bowtie2 -x genome_index \
  -U trimmed_sample.fq.gz \
  --very-sensitive -k 1 \
  -S sample_aligned.sam 2> sample_align.log
```

**植物比对注意事项**:
- 推荐 `--very-sensitive` 模式提高比对率
- `-k 1` 限制唯一比对，降低重复序列干扰
- 对多倍体物种，考虑使用亚基因组特异性比对

### 1.3 SAM 到 BAM 转换与排序

```bash
samtools view -bS -q 30 sample_aligned.sam > sample_quality.bam
samtools sort -o sample_sorted.bam sample_quality.bam
samtools index sample_sorted.bam
```

**MAPQ 过滤**: 推荐 MAPQ >= 30 (唯一比对)

### 1.4 去 PCR 重复

```bash
picard MarkDuplicates \
  I=sample_sorted.bam \
  O=sample_dedup.bam \
  M=sample_dedup_metrics.txt \
  REMOVE_DUPLICATES=true

samtools index sample_dedup.bam
```

### 1.5 ChIP-seq 特殊质控

```bash
# Fingerprint plot 评估 ChIP 信号富集度
plotFingerprint -b ChIP_dedup.bam Input_dedup.bam \
  --labels ChIP Input \
  -o qc/fingerprint.pdf

# 评估片段长度分布
samtools view ChIP_dedup.bam | \
  awk '{print length($10)}' | sort -n | uniq -c \
  > qc/fragment_length_distribution.txt
```

**质量指标**:
- **唯一比对率**: > 70%
- **Fingerprint enrichment**: ChIP 曲线明显偏离 Input 对角线
- **FRiP (Fraction of Reads in Peaks)**: 转录因子 > 5%, 组蛋白修饰 > 20%

---

## 2. Peak Calling

### 2.1 转录因子 ChIP-seq (Narrow Peak)

```bash
macs2 callpeak \
  -t TF_ChIP_dedup.bam \
  -c Input_dedup.bam \
  -f BAM -g 1.19e8 \
  -n TF_sample \
  --nomodel --extsize 200 \
  --outdir macs2_output/TF_sample/

# 查看 peak 数量
wc -l macs2_output/TF_sample/TF_sample_peaks.narrowPeak
```

**拟南芥基因组大小**: `-g 1.19e8`
**水稻基因组大小**: `-g 3.74e8`
**玉米基因组大小**: `-g 2.04e9`
**大豆基因组大小**: `-g 9.75e8`

### 2.2 组蛋白修饰 ChIP-seq (Broad Peak)

```bash
macs2 callpeak \
  -t H3K27ac_ChIP_dedup.bam \
  -c Input_dedup.bam \
  -f BAM -g 1.19e8 \
  -n H3K27ac_sample \
  --broad --broad-cutoff 0.1 \
  --outdir macs2_output/H3K27ac_sample/
```

**常见组蛋白修饰及 peak 类型**:

| 修饰 | Peak 类型 | 生物学意义 |
|------|----------|-----------|
| H3K4me3 | Narrow | 活跃启动子 |
| H3K27ac | Broad | 活跃增强子/启动子 |
| H3K27me3 | Broad | 抑制性染色质 |
| H3K9me2 | Broad | 转座子沉默 (植物特有) |
| H3K36me3 | Broad | 转录延伸 |

---

## 3. Peak 注释

### 3.1 ChIPseeker 注释

```r
library(ChIPseeker)
library(clusterProfiler)

# 载入基因组注释 (拟南芥示例)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

# 对于其他植物物种，使用 GFF 文件
# txdb <- makeTxDbFromGFF("Oryza_sativa.gff3")

# 读取 peaks
peaks <- readPeakFile("tf_peaks.narrowPeak")

# 注释 peaks
peak_anno <- annotatePeak(peaks,
                          TxDb = txdb,
                          tssRegion = c(-3000, 3000),
                          annoDb = "org.At.tair.db")

# 可视化
plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno)

# 导出注释结果
write.csv(as.data.frame(peak_anno), "peak_annotation.csv")
```

### 3.2 靶基因 GO 富集分析

```r
# 提取靶基因
target_genes <- unique(as.data.frame(peak_anno)$geneId)

# GO 富集
library(clusterProfiler)
ego <- enrichGO(gene = target_genes,
                OrgDb = org.At.tair.db,
                keyType = "TAIR",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# 可视化
dotplot(ego)
barplot(ego)
```

---

## 4. Motif 富集分析

### 4.1 提取 peak 中心序列

```bash
# 提取 peak summit 周围 +/- 100bp 序列
bedtools slop -i peaks_summits.bed -g chrom.sizes -b 100 \
  > summit_regions.bed
bedtools getfasta -fi genome.fa -bed summit_regions.bed \
  -fo summit_sequences.fa
```

### 4.2 MEME-ChIP Motif 发现

```bash
meme-chip -oc meme_output/ \
  -db JASPAR2022_CORE_plants_non_redundant.meme \
  -meme-mod anr -nmeme 5 \
  summit_sequences.fa
```

---

## 5. 差异结合分析

### 5.1 DiffBind 分析

```r
library(DiffBind)

# 读取 samplesheet
samples <- read.csv("samplesheet.csv")
dba_obj <- dba(sampleSheet = samples)

# 计数
dba_obj <- dba.count(dba_obj, bUseSummarizeOverlaps = TRUE)

# 差异分析
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION,
                         minMembers = 2)
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

# 结果
results <- dba.report(dba_obj, method = DBA_DESEQ2)

# 可视化
dba.plotMA(dba_obj)
dba.plotPCA(dba_obj, attributes = DBA_CONDITION)
dba.plotVolcano(dba_obj)
```

---

## 6. 可视化

### 6.1 基因组浏览器 Track

```bash
# 生成 bigWig 用于 IGV
bamCoverage -b sample_dedup.bam -o sample_RPKM.bw \
  --normalizeUsing RPKM --binSize 10
```

### 6.2 TSS 附近信号分布

```bash
computeMatrix reference-point -S sample_RPKM.bw \
  -R genes.bed --referencePoint TSS \
  -a 3000 -b 3000 -o matrix_TSS.gz

plotProfile -m matrix_TSS.gz -o profile_TSS.pdf --perGroup
plotHeatmap -m matrix_TSS.gz -o heatmap_TSS.pdf --colorMap Blues
```

---

## 常见问题

### Q: 比对率低怎么办?
A: 植物基因组通常较大且多倍体化，尝试使用 `--very-sensitive` 模式。检查是否去除了接头序列。考虑允许 multi-mapping (`-k 2` 或 `-k 3`) 用于重复序列丰富的区域。

### Q: Peak 数量太少?
A: 降低 MACS2 的 q-value 阈值 (默认 0.05，可尝试 0.1)。检查 IP 效率和比对质量。考虑排除黑名单区域（如果有）。

### Q: Input 对照信号异常高?
A: 检查 Input DNA 浓度和超声打断效率。考虑使用 IgG 对照替代或补充。

### Q: 植物特有组蛋白修饰的分析策略?
A: H3K9me2 和 H3K27me3 在植物中分布广泛，需要使用 broad peak 模式。注意排除转座子/重复区域的假阳性信号。

---

## 参考

- MACS2: https://github.com/macs3-project/MACS
- deepTools: https://deeptools.readthedocs.io/
- ChIPseeker: https://bioconductor.org/packages/ChIPseeker/
- DiffBind: https://bioconductor.org/packages/DiffBind/
- MEME-ChIP: https://meme-suite.org/meme/tools/meme-chip
- HOMER: http://homer.ucsd.edu/homer/
