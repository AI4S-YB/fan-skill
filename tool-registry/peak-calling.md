# ATAC-seq Peak Calling -- 工具目录

## 概述

Peak calling 是 ATAC-seq 分析的核心步骤，用于识别开放染色质区域。ATAC-seq 的数据特征与 ChIP-seq 不同，需要专门的参数设置。

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| qvalue (FDR) | 0.05 | Low signal-to-noise data: relax to 0.10; high-quality replicates with strong enrichment: tighten to 0.01 | ATAC-seq in plants can have higher background due to organelle DNA contamination; strict cutoff may miss genuine open regions in low-coverage datasets |
| extsize | 200 | Nucleosome-free regions focus: use 100; nucleosome positioning analysis: use 147 | Tn5 cuts ~9 bp from binding site; 200 bp captures the protected region between two adjacent Tn5 insertion events |
| shift | -100 | Paired-end data with actual fragment lengths: set to 0; Tn5 calibration experiment available: adjust empirically | Compensates for the known ~4-5 bp offset of Tn5 insertion on each strand of the cut site relative to the binding midpoint |
| --nomodel | always on | Never turn off | Tn5 transposase cutting is fundamentally different from sonication; MACS2's model assumes ChIP-seq fragment distribution which does not apply to ATAC-seq |
| --keep-dup | all | Very deep sequencing (>100M mapped reads): consider "auto" | ATAC-seq open chromatin regions legitimately produce reads at identical positions due to Tn5 insertion preference; removing duplicates loses signal in highly accessible regions |

## 推荐工具

### 1. MACS2 (ATAC-seq 模式) -- 金标准

**描述**: MACS2 是 ATAC-seq peak calling 的首选工具。使用 `--nomodel --shift -100 --extsize 200` 参数组合补偿 Tn5 转座酶的切割偏好。

**适用场景**: 所有植物 ATAC-seq 数据的 peak calling。

**关键参数**:
- `-t <bam>`: 处理样本 BAM 文件（ATAC-seq 无 control，仅需单样本）
- `-f BAM`: 输入文件格式
- `-g <genome_size>`: 有效基因组大小
  - 拟南芥: `-g 1.19e8`
  - 水稻: `-g 3.74e8`
  - 玉米: `-g 2.04e9`
  - 大豆: `-g 9.75e8`
  - 番茄: `-g 7.81e8`
- `-n <name>`: 输出前缀
- `--nomodel`: 必须使用 -- Tn5 酶切割模式不同于超声打断
- `--shift -100`: Tn5 结合位点偏移校正
- `--extsize 200`: 延伸 200bp 模拟核小体间区域
- `--keep-dup all`: 保留所有可能的重复 reads
- `-q 0.05`: FDR 阈值
- `--outdir <dir>`: 输出目录

**基本示例**:
```bash
macs2 callpeak -t sample_dedup.bam \
  -f BAM -g 1.19e8 \
  -n sample_atac \
  --nomodel --shift -100 --extsize 200 \
  --keep-dup all -q 0.05 \
  --outdir macs2_output/
```

**不加 Input 对照**:
ATAC-seq 通常不需要 Input 对照，因为 Tn5 切割产生的是开放染色质信号而非特定蛋白质富集信号。

**输出文件**:
- `*_peaks.narrowPeak`: Peak 位置 BED 文件
- `*_peaks.xls`: Peak 统计表格
- `*_summits.bed`: Peak 峰顶位置（用于 motif 分析）
- `*.bedGraph`: 信号覆盖文件

### 2. Genrich (备选)

**描述**: Genrich 专为 ATAC-seq peak calling 设计，可处理多重复数据。

**基本用法**:
```bash
Genrich -t sample1.bam,sample2.bam,sample3.bam \
  -o genrich_peaks.narrowPeak \
  -f genrich_log.txt \
  -j -y -r -v
```

**参数说明**:
- `-t`: 多个处理 BAM 文件（逗号分隔）
- `-j`: ATAC-seq 模式
- `-y`: 保留重复 reads
- `-r`: 移除 PCR 重复
- `-v`: 详细输出

---

## 植物 ATAC-seq Peak Calling 注意事项

1. **不使用 `--broad` 模式**: ATAC-seq peaks 是窄而尖的信号峰，不应使用 broad peak 模式
2. **基因组大小**: 必须使用有效基因组大小（`mappable genome size`），不是总基因组大小
3. **去重复**: 推荐保留所有 reads（`--keep-dup all`），因为 ATAC-seq 中相同位置的 reads 可能是真信号
4. **线粒体/叶绿体**: peak calling 前应移除线粒体和叶绿体 reads
5. **覆盖度**: 对植物基因组，推荐每个样本至少 25M uniquely mapped reads
6. **预期 peak 数**:
   - 拟南芥叶片: ~20,000-40,000 peaks
   - 水稻叶片: ~30,000-60,000 peaks
   - 玉米叶片: ~50,000-100,000 peaks

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Too few peaks called (<5,000) | Stringent qvalue or low sequencing coverage | Relax -q to 0.10; verify >=25M uniquely mapped reads per sample after organelle removal |
| Peaks enriched in organelle genome coordinates | Mitochondria/chloroplast reads not removed before peak calling | Filter BAM to nuclear chromosomes only before peak calling; use `samtools view -L nuclear_regions.bed` |
| MACS2 crashes with "pair-ended" error | BAM not coordinate-sorted or index missing | `samtools sort -o sorted.bam input.bam && samtools index sorted.bam` |
| Broad, flat signal without discrete narrow peaks | Low Tn5 efficiency or over-fixed tissue | Check Tn5 enzyme activity with positive control; verify nuclei isolation protocol; consider adding 0.1% digitonin |
