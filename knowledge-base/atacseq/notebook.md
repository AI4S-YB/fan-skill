# ATAC-seq 数据分析 -- 分析笔记本

## 分析概览

本笔记本提供植物 ATAC-seq 数据的完整分析流程，涵盖开放染色质区域的鉴定、差异可及性分析、转录因子足迹分析和功能注释。

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
- GC content: 与参考基因组 GC 含量接近
- Adapter content: ATAC-seq 使用 Nextera 转座酶接头，需检查 Tn5 接头污染

### 1.2 接头与质量修剪

```bash
# Trimmomatic 去接头
trimmomatic PE -threads 8 \
  sample_R1.fq.gz sample_R2.fq.gz \
  sample_R1_trimmed.fq.gz sample_R1_unpaired.fq.gz \
  sample_R2_trimmed.fq.gz sample_R2_unpaired.fq.gz \
  ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:keepBothReads \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 1.3 Reads 比对到参考基因组

```bash
# Bowtie2 比对 (植物推荐 --very-sensitive)
bowtie2 -x genome_index \
  -1 sample_R1_trimmed.fq.gz -2 sample_R2_trimmed.fq.gz \
  --very-sensitive -X 2000 --no-discordant --no-mixed \
  -S sample_aligned.sam 2> sample_align.log
```

**植物比对注意事项**:
- 推荐 `--very-sensitive` 模式提高比对率
- `-X 2000` 允许更大插入片段（植物核小体间距约 200bp）
- 对多倍体物种，考虑亚基因组特异性比对

### 1.4 SAM/BAM 处理

```bash
# 过滤低质量比对
samtools view -bS -q 30 -F 4 sample_aligned.sam > sample_filtered.bam

# 排序和索引
samtools sort -o sample_sorted.bam sample_filtered.bam
samtools index sample_sorted.bam
```

### 1.5 去线粒体和叶绿体 reads

```bash
# 提取核基因组 reads (拟南芥示例)
samtools view -b sample_sorted.bam chr1 chr2 chr3 chr4 chr5 > sample_nuclear.bam
samtools index sample_nuclear.bam
```

**植物特有**: 绿色组织 ATAC-seq 中叶绿体 reads 比例可能很高（>50%），需去除。

### 1.6 去 PCR 重复

```bash
picard MarkDuplicates \
  I=sample_nuclear.bam \
  O=sample_dedup.bam \
  M=sample_dedup_metrics.txt \
  REMOVE_DUPLICATES=true

samtools index sample_dedup.bam
```

### 1.7 ATAC-seq 特有质控

```bash
# 插入片段长度分布 - ATAC-seq 关键质控
samtools view sample_dedup.bam | \
  awk '{print sqrt($9^2)}' | sort -n | uniq -c \
  > qc/fragment_length_distribution.txt

# 计算 TSS 富集分数
computeMatrix reference-point -S sample_RPKM.bw \
  -R genes.bed --referencePoint TSS \
  -a 2000 -b 2000 -o matrix_TSS.gz

plotProfile -m matrix_TSS.gz -o profile_TSS.pdf
```

**ATAC-seq 特有质量指标**:
- **Nucleosome-free region (< 150bp) 比例**: 应有明显峰
- **Mono-nucleosome (~180-247bp)**: 弱峰
- **Di-nucleosome (~315-473bp)**: 更弱
- **TSS 富集分数 (TSS enrichment)**: > 7 为良好
- **FRiP (Fraction of Reads in Peaks)**: > 30%

---

## 2. Peak Calling

### 2.1 MACS2 ATAC-seq 模式

ATAC-seq 的数据特征不同于 ChIP-seq：转座酶切割产生的是开放染色质区域信号，而非蛋白质结合峰。

```bash
macs2 callpeak \
  -t sample_dedup.bam \
  -f BAM -g 1.19e8 \
  -n sample_atac \
  --nomodel --shift -100 --extsize 200 \
  --keep-dup all \
  --outdir macs2_output/

# 查看 peak 数量
wc -l macs2_output/sample_atac_peaks.narrowPeak
```

**ATAC-seq Peak Calling 关键参数**:
- `--nomodel`: ATAC-seq 必须使用此参数（Tn5 切割模式不同于超声打断）
- `--shift -100`: Tn5 切割偏移校正
- `--extsize 200`: 延伸 200bp（模拟核小体间距）
- `--keep-dup all`: 保留所有 reads（ATAC-seq 不推荐去重）
- `-q 0.05`: FDR 阈值

**植物基因组大小参考**:
| 物种 | `-g` 参数 |
|------|----------|
| 拟南芥 | 1.19e8 |
| 水稻 | 3.74e8 |
| 玉米 | 2.04e9 |
| 大豆 | 9.75e8 |

---

## 3. Peak 注释

### 3.1 ChIPseeker 注释

```r
library(ChIPseeker)
library(clusterProfiler)

# 对于拟南芥
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

# 对于其他植物，使用 GFF 文件
# txdb <- makeTxDbFromGFF("Oryza_sativa.gff3")

# 读取 peaks
peaks <- readPeakFile("sample_atac_peaks.narrowPeak")

# 注释 peaks
peak_anno <- annotatePeak(peaks,
                          TxDb = txdb,
                          tssRegion = c(-3000, 3000),
                          annoDb = "org.At.tair.db")

# 可视化
plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno)

write.csv(as.data.frame(peak_anno), "peak_annotation.csv")
```

### 3.2 基因组区域分布

ATAC-seq peaks 在基因组中的区域分布反映了开放染色质的调控特征：

- **Promoter (<=1kb)**: 与活跃转录直接相关
- **Promoter (1-2kb)**: 远端调控
- **5' UTR**: 转录调控
- **3' UTR**: 转录后调控
- **Exon / Intron**: 可能涉及基因调控
- **Intergenic**: 增强子或非编码调控元件

---

## 4. 差异可及性分析

### 4.1 DiffBind + DESeq2

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

# 导出
write.csv(as.data.frame(results), "differential_accessibility.csv")
```

**重要提示**:
- 需要至少 2 个生物学重复
- 使用 `bUseSummarizeOverlaps=TRUE` 更准确的计数
- 推荐 DESeq2 方法（适用于低重复情况）

---

## 5. TF Footprinting

### 5.1 TOBIAS Footprinting

```bash
# 步骤 1: ATACorrect (偏差校正)
TOBIAS ATACorrect --bam sample_dedup.bam \
  --genome genome.fa --peaks sample_peaks.narrowPeak \
  --outdir tobias_output/

# 步骤 2: ScoreMotifs (motif footprint 评分)
TOBIAS ScoreMotifs --signals tobias_output/sample_corrected.bw \
  --motifs JASPAR2022_CORE_plants_non_redundant.meme \
  --genome genome.fa --outdir tobias_output/motif_scores/

# 步骤 3: BINDetect (差异 footprint 检测)
TOBIAS BINDetect \
  --signals condition1_corrected.bw condition2_corrected.bw \
  --motifs JASPAR2022_CORE_plants_non_redundant.meme \
  --genome genome.fa \
  --conditions condition1 condition2 \
  --outdir tobias_output/differential/
```

**植物特有 TF 数据库**:
- PlantTFDB: http://planttfdb.gao-lab.org/
- JASPAR Plants: https://jaspar.genereg.net/
- CIS-BP: http://cisbp.ccbr.utoronto.ca/

---

## 6. Motif 分析

### 6.1 HOMER Motif 富集

```bash
# 提取 peak 中心序列
bedtools slop -i sample_peaks.narrowPeak -g chrom.sizes -b 100 \
  > summit_regions.bed
bedtools getfasta -fi genome.fa -bed summit_regions.bed \
  -fo summit_sequences.fa

# HOMER motif 富集分析
findMotifsGenome.pl summit_regions.bed genome.fa homer_output/ \
  -size given -mask
```

---

## 7. 可视化

### 7.1 基因组浏览器 Track

```bash
# 生成 bigWig
bamCoverage -b sample_dedup.bam -o sample_CPM.bw \
  --normalizeUsing CPM --binSize 10
```

### 7.2 可及性热图

```bash
# TSS 区域可及性热图
computeMatrix reference-point -S sample_CPM.bw \
  -R genes.bed --referencePoint TSS \
  -a 2000 -b 2000 -o matrix_TSS.gz

plotHeatmap -m matrix_TSS.gz -o heatmap_TSS.pdf \
  --colorMap YlOrRd --zMax 10
```

---

## 常见问题

### Q: ATAC-seq 比对率偏低怎么办？
A: 检查线粒体/叶绿体 reads 比例。植物绿色组织中 cpDNA reads 可达 50%+，此为正常现象，去除后再计算比对率。尝试 `--very-sensitive-local` 模式。

### Q: Peak 数目异常多或异常少？
A: 正常植物 ATAC-seq peak 数因物种和组织而异：拟南芥叶片 ~20,000-40,000；水稻叶片 ~30,000-60,000；玉米叶片 ~50,000-100,000。

### Q: Fragment length distribution 无核小体 ladder？
A: 检查 Tn5 酶切效率。可能需要调整细胞核分离和转座反应条件。也可以放宽 MACS2 参数。

### Q: 差异可及性分析时无结果？
A: 检查生物学重复的一致性（PCA 图）。考虑使用 Consensus peaks（至少 2 个重复共有）而非单个重复的 peaks。

---

## 植物 ATAC-seq 实验的特有挑战

### 细胞壁去除与核分离

植物 ATAC-seq 最关键也最容易出错的步骤是核分离：

**为什么细胞壁是问题**：
- Tn5 转座酶无法穿透细胞壁，必须先去除细胞壁释放细胞核
- 细胞壁去除不充分 → 核释放不足 → Tn5 富集信号弱
- 细胞壁去除过度 → 细胞器膜破裂 → 叶绿体和线粒体 DNA 污染

**常用核分离方法对比**：

| 方法 | 适用组织 | 优势 | 劣势 |
|------|---------|------|------|
| 刀片切碎 + 过滤 | 叶片、幼苗 | 简单、快速、无需特殊设备 | 对纤维组织效果差 |
| 液氮研磨 + 蔗糖梯度 | 根、茎、种子 | 核纯度高 | 耗时、需液氮和超速离心 |
| 渗透裂解 + NP-40 | 原生质体 | 核完整性好 | 需先制备原生质体 |
| 压力破碎（French press） | 木质化组织 | 对硬组织有效 | 需要专门设备 |
| 匀浆器（Dounce homogenizer） | 软嫩组织（花、幼叶） | 温和、适合少量材料 | 处理量有限 |

**关键 QC 指标**：
- 核数量和完整性：台盼蓝或 DAPI 染色后在显微镜下计数。需 50,000-100,000 个完整核用于标准 ATAC-seq
- 核纯度：无可见细胞碎片和叶绿体（叶绿体在 DAPI 下呈红色自发光）
- 裂解效率：> 80% 细胞裂解是目标

### 不同植物组织的特殊处理

**叶片**:
- 最常见但叶绿体污染严重（50-80% cpDNA reads）
- 建议使用幼嫩叶片（叶绿体较小较少）而非成熟叶片
- 刀片切碎法通常效果最好（对叶肉细胞）
- 如果 cpDNA reads 持续过高，考虑在缓冲液中加入 0.5% Triton X-100 帮助去除叶绿体膜

**根**:
- 根组织富含纤维，刀片切碎效果差，液氮研磨更合适
- 根的内皮层（凯氏带）可能阻碍缓冲液渗透
- 根表面微生物污染需要先清洗去除
- 根尖分生组织细胞核密度高，是 ATAC-seq 的良好材料

**花/花序**:
- 花器官细胞类型多样，ATAC 信号是多种细胞类型的平均值
- 花粉粒含有高度压缩的精细胞核，其 ATAC 信号特征不同于体细胞，可能引入异质性
- 对于禾本科（水稻、玉米、小麦），花序结构复杂，建议解剖后选取特定部位（如花药或子房）

**种子/胚**:
- 种子含有大量储存蛋白和油脂，干扰核提取
- 吸胀（imbibition）后的种子核提取效率更高
- 成熟干燥种子几乎无法进行 ATAC-seq，需使用发育中的胚

### 田间 vs 温室样本

植物 ATAC-seq 中环境条件对染色质可及性有显著影响：

**温室样本**（推荐用于首次 ATAC-seq）：
- 条件可控（光照、温度、湿度）
- 生物学重复间变异性小
- 适合研究特定处理效应（如激素处理、胁迫）
- 取样时间和光照条件可控

**田间样本**：
- 条件不可控但更具生态真实性
- 生物学重复间变异性可能很大（微环境差异）
- 建议增加重复数（n >= 5，而非温室的 n >= 3）
- 必须精确记录取样时的环境参数（光照、温度、土壤湿度、生长阶段）
- 田间植物的染色质状态受昼夜节律、温度波动、UV 辐射等多重因素影响
- 单个田间时间点可能是特殊天气条件（如取样前一日高温）的"快照"，不代表该生长阶段的典型状态

### 非模式植物的特异性考虑

大多数植物 ATAC-seq 目前仍集中在拟南芥和水稻上。对于非模式物种：

- **无参考基因组**：如果近缘物种有高质量参考（同属或同科，<20Mya 分化），可尝试 cross-species mapping。比对率可能降低但 ATAC 信号富集区域仍是可解释的。使用较宽松的比对参数（--very-sensitive-local）提高 cross-species 比对率。
- **倍性**：多倍体物种中 peak calling 需要考虑 subgenome specificity。如果已有亚基因组组装，分别比对。如果全基因组未分型，peak calling 后注意保守 peak 的 counts 可能来源自多个亚基因组。
- **GC 含量**：许多植物基因组 GC 含量偏低（30-40%），而 Tn5 有轻微的 GC 偏好。在 footprinting 分析中使用 TOBIAS 的 ATACorrect 步骤校正此类偏差。

---

## 参考

- MACS2: https://github.com/macs3-project/MACS
- deepTools: https://deeptools.readthedocs.io/
- ChIPseeker: https://bioconductor.org/packages/ChIPseeker/
- DiffBind: https://bioconductor.org/packages/DiffBind/
- TOBIAS: https://github.com/loosolab/TOBIAS
- HOMER: http://homer.ucsd.edu/homer/
- ATAC-seq Guidelines: https://www.encodeproject.org/atac-seq/
