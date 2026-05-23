# ChIP-seq 数据预处理 — 工具目录

## 概述

ChIP-seq 数据预处理包括原始数据质量控制、reads 比对到参考基因组、去重复和格式转换。由于植物基因组的特殊性（大基因组、多倍体、重复序列），预处理需要特别注意。

## 推荐工具

### 1. FastQC (质量评估)

**描述**: 测序数据质量评估的标准工具。

**使用**:
```bash
fastqc sample_chip.fastq.gz -o qc_output/
multiqc qc_output/ -o qc_report/
```

### 2. Bowtie2 (比对，推荐)

**描述**: 适用于短 reads 比对的常用工具，ChIP-seq 分析的首选比对工具。

**关键参数**:
- `-x <index>`: 参考基因组索引
- `-U <reads>`: 单端测序 reads
- `-1/-2 <reads>`: 双端测序 reads
- `--very-sensitive`: 植物基因组推荐使用灵敏模式
- `-k 1`: 仅报告唯一比对 (植物基因组重复序列较多)

**示例**:
```bash
bowtie2 -x genome_index -U trimmed_sample.fq.gz \
  --very-sensitive -k 1 \
  -S sample_aligned.sam 2> alignment_stats.log
```

### 3. BWA-MEM (备选)

**描述**: 适用于较长 reads (>= 70bp) 的比对工具。

**适用场景**: 双端测序 reads > 100bp 时推荐使用 BWA-MEM。

### 4. SAMtools (SAM/BAM 处理)

**描述**: SAM/BAM 文件格式转换、排序、索引。

**预处理流程**:
```bash
# SAM -> BAM
samtools view -bS sample_aligned.sam > sample.bam

# 按坐标排序
samtools sort -o sample_sorted.bam sample.bam

# 建立索引
samtools index sample_sorted.bam
```

### 5. Picard / Sambamba (去重复)

**描述**: 去除 PCR 重复 reads，减少 PCR 扩增偏差。

```bash
picard MarkDuplicates \
  I=sample_sorted.bam \
  O=sample_dedup.bam \
  M=dedup_metrics.txt \
  REMOVE_DUPLICATES=true
```

### 6. deepTools (质量控制)

**描述**: 用于 ChIP-seq 数据质量评估的综合工具集。

**关键分析**:
- `plotFingerprint`: 评估 ChIP 信号富集度
- `plotCoverage`: 基因组覆盖度
- `multiBamSummary`: 多样本相关性分析

### 植物特殊注意事项

- 使用 `--very-sensitive` 模式提高比对率
- 使用 `-k 1` 限制唯一比对，避免重复序列干扰
- 大基因组物种 (玉米、小麦) 注意内存使用
- 多倍体物种考虑使用亚基因组特异比对
