# ATAC-seq 数据预处理 -- 工具目录

## 概述

ATAC-seq 数据预处理包括原始数据质控、接头修剪、基因组比对、低质量过滤和去重复。植物数据还需要特别处理线粒体和叶绿体来源的 reads。

## 推荐流程

### 步骤 1: 质控 (FastQC + MultiQC)

```bash
# 对原始数据进行质量评估
fastqc -t 8 raw/*.fastq.gz -o qc/raw/
multiqc qc/raw/ -o qc/raw_report/
```

### 步骤 2: 接头与质量修剪 (Trimmomatic / Cutadapt)

ATAC-seq 的 Nextera Tn5 接头序列:
- Read 1: CTGTCTCTTATACACATCT
- Read 2: CTGTCTCTTATACACATCT

```bash
# 使用 Trimmomatic
trimmomatic PE -threads 8 \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  sample_R1_paired.fastq.gz sample_R1_unpaired.fastq.gz \
  sample_R2_paired.fastq.gz sample_R2_unpaired.fastq.gz \
  ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:keepBothReads \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# 或使用 Cutadapt (更快)
cutadapt -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -o sample_R1_trimmed.fastq.gz \
  -p sample_R2_trimmed.fastq.gz \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  -m 36 -j 8
```

**参数说明**:
- `LEADING:3 / TRAILING:3`: 切除首尾低质量碱基（Q<3）
- `SLIDINGWINDOW:4:15`: 滑动窗口平均质量 < 15 时截断
- `MINLEN:36`: 丢弃 < 36bp 的 reads

### 步骤 3: 基因组比对 (Bowtie2 / BWA)

**Bowtie2 (推荐)**:
```bash
# 构建基因组索引
bowtie2-build genome.fa genome_index

# 比对
bowtie2 -x genome_index \
  -1 sample_R1_trimmed.fastq.gz \
  -2 sample_R2_trimmed.fastq.gz \
  --very-sensitive \
  -X 2000 --no-discordant --no-mixed \
  -p 16 \
  -S sample_aligned.sam 2> sample_align.log
```

**BWA-MEM (备选，适合大基因组)**:
```bash
bwa mem -t 16 genome.fa \
  sample_R1_trimmed.fastq.gz sample_R2_trimmed.fastq.gz \
  > sample_aligned.sam
```

**比对参数说明**:
| 参数 | 说明 |
|------|------|
| `--very-sensitive` | 最大灵敏度模式，提高比对率 |
| `-X 2000` | 最大插入片段长度（植物核小体间距约 200bp） |
| `--no-discordant` | 只保留合理配对的 reads |
| `--no-mixed` | 不允许单端比对 |

### 步骤 4: BAM 后处理

```bash
# 过滤 MAPQ < 30 的低质量比对
samtools view -bS -q 30 -F 4 sample_aligned.sam > sample_filtered.bam

# 排序
samtools sort -@ 8 -o sample_sorted.bam sample_filtered.bam

# 索引
samtools index sample_sorted.bam
```

### 步骤 5: 去除线粒体和叶绿体 reads (植物特有)

```bash
# 方法 1: 按染色体名排除
# 拟南芥: 线粒体=chrM, 叶绿体=chrC
samtools view -b sample_sorted.bam \
  chr1 chr2 chr3 chr4 chr5 > sample_nuclear.bam

# 方法 2: 使用 SAMtools region
# 先获取核基因组 region
grep ">" genome.fa | sed 's/>//' | grep -v -E "Mt|Pt|ChrC|ChrM" \
  > nuclear_chromosomes.txt
samtools view -b sample_sorted.bam \
  -L nuclear_chromosomes.txt > sample_nuclear.bam

# 构建索引
samtools index sample_nuclear.bam
```

### 步骤 6: 去 PCR 重复

```bash
# Picard MarkDuplicates
picard MarkDuplicates \
  I=sample_nuclear.bam \
  O=sample_dedup.bam \
  M=sample_dedup_metrics.txt \
  REMOVE_DUPLICATES=true \
  VALIDATION_STRINGENCY=LENIENT

# 或使用 SAMtools (更快)
samtools markdup -r sample_nuclear.bam sample_dedup.bam
samtools index sample_dedup.bam
```

### 步骤 7: 片段大小校正 (Tn5 偏移)

```bash
# deepTools alignmentSieve 校正 Tn5 偏移
alignmentSieve \
  --bam sample_dedup.bam \
  --filterMetrics fragment_length_metrics.txt \
  --ATACshift \
  -o sample_shifted.bam
```

### 步骤 8: 最终质控

```bash
# 统计比对结果
samtools flagstat sample_dedup.bam > sample_flagstat.txt

# 插入片段长度分布 (ATAC-seq 关键指标)
samtools view sample_dedup.bam | \
  awk '{print sqrt($9^2)}' | sort -n | uniq -c \
  > fragment_sizes.txt

# 生成 bigWig 用于可视化
bamCoverage -b sample_dedup.bam \
  --normalizeUsing CPM \
  --binSize 10 \
  -o sample_CPM.bw
```

---

## 植物特有考虑

### 线粒体和叶绿体 reads

植物 ATAC-seq 中最显著的污染源。绿色组织中 cpDNA 可达总 reads 的 50%+，根部组织中 mtDNA 比例较高。

**预期比例**:
| 组织 | 核 DNA | 叶绿体 DNA | 线粒体 DNA |
|------|--------|-----------|-----------|
| 绿色叶片 | 40-60% | 30-50% | 5-10% |
| 根 | 70-80% | <5% | 15-25% |
| 花 | 60-80% | 10-20% | 10-15% |
| 愈伤组织 | 80-95% | <5% | <5% |

### 多倍体比对策略

对于多倍体植物（如小麦 6x），比对时将有挑战：
- 使用 `-k 1` 限制唯一比对可能丢失同源区域信号
- 使用 `-k 3` 允许部分 multi-mapping
- 考虑亚基因组特异性参考基因组

---

## 质量指标

| 指标 | 优秀 | 可接受 | 需要关注 |
|------|------|--------|---------|
| 总 reads 数 | >50M | 25-50M | <25M |
| 比对率 | >90% | 70-90% | <70% |
| 核 DNA 比例 | >80% | 50-80% | <50% |
| PCR Duplicate % | <20% | 20-40% | >40% |
| Fragment NFR % | >60% | 40-60% | <40% |
