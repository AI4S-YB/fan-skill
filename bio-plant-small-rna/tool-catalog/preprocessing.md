# 小 RNA 数据预处理 — 工具目录

## 概述

植物小 RNA 测序数据预处理包括质量评估、接头去除、长度筛选和低质量 reads 过滤。

## 推荐工具

### 1. FastQC (质量评估)

**描述**: 测序数据质量评估的标准工具。

**使用**:
```bash
fastqc sample.fastq.gz -o qc_output/
```

### 2. Cutadapt (接头去除)

**描述**: 去除测序接头序列。

**关键参数**:
- `-a <adapter>`: 3' 端接头序列
- `-m <min_len>`: 最小 reads 保留长度 (推荐 18)
- `-M <max_len>`: 最大 reads 保留长度 (推荐 30)
- `--discard-untrimmed`: 丢弃未检测到接头的 reads

**植物小 RNA 常见接头**: `TGGAATTCTCGGGTGCCAAGG`

**示例**:
```bash
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 30 \
  -o trimmed.fastq.gz raw.fastq.gz
```

### 3. FASTX-Toolkit (长度筛选和质控)

**描述**: 用于 FASTQ/A 文件处理的工具集。

**关键工具**:
- `fastq_quality_filter`: 质量过滤
- `fastq_to_fasta`: 格式转换

### 预处理流程

```
Raw FASTQ → FastQC → Cutadapt (去接头) → 长度筛选 (18-30nt) → 
质量过滤 (Q20) → 去除 rRNA/tRNA → Clean reads
```
