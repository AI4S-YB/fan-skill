# 组装后碱基校正 (Polishing) -- 工具目录

## 概述

组装后碱基校正（Polishing）用于修复组装中残留的序列错误。不同测序平台需要不同的 polishing 工具：ONT 数据用 Medaka，PacBio HiFi 数据用 gcpp/Arrow，Illumina 短读长数据用 Pilon。

## 推荐工具

### 1. Medaka -- ONT Polishing 专用

**描述**: Medaka 是 ONT 开发的神经网络碱基校正工具，利用原始信号或 reads 比对来修正组装错误。

**安装**:
```bash
conda install -c bioconda medaka
```

**标准流程**:
```bash
# 步骤 1: 将 ONT reads 比对到组装
mini_align -i ont_reads.fastq.gz \
  -r assembly.fasta \
  -m -p alignment -t 16

# 步骤 2: Medaka consensus
medaka_consensus -i ont_reads.fastq.gz \
  -d assembly.fasta \
  -o medaka_output/ \
  -t 16 \
  -m r941_min_hac_g507
```

**模型选择** (`-m` 参数):
| 模型 | 适用数据 |
|------|---------|
| `r941_min_hac_g507` | ONT MinION R9.4.1, HAC 模式 |
| `r941_min_sup_g507` | ONT MinION R9.4.1, SUP 模式 |
| `r1041_e82_260bps_hac_v4.1.0` | ONT R10.4.1, 260 bps |
| `r1041_e82_400bps_hac_v4.1.0` | ONT R10.4.1, 400 bps |
| `r10_min_high_g360` | ONT MinION R10 |

**运行时间和内存**:
- 内存: 8-16 GB (较小基因组) 到 128 GB+ (大基因组)
- 运行时间: 与基因组大小线性相关

### 2. gcpp (Genomic Consensus) -- PacBio HiFi Polishing

**描述**: gcpp (原 Arrow) 是 PacBio 的 HiFi 数据 polishing 工具，利用 CCS reads 的多序列比对来提高碱基准确性。

**安装**:
```bash
conda install -c bioconda pbgcpp
```

```bash
# 步骤 1: 将 HiFi reads 比对到组装
pbmm2 align hifi_reads.bam assembly.fasta \
  hifi_aligned.bam --sort -j 16

# 步骤 2: gcpp polishing
gcpp -j 16 \
  -r assembly.fasta \
  -o polished.fasta \
  hifi_aligned.bam
```

### 3. Pilon -- Illumina 短读长 Polishing

**描述**: Pilon 使用 Illumina 短读长的高碱基准确性（Q30+）来校正组装中的小规模错误（SNP 和短 indel）。

**安装**:
```bash
conda install -c bioconda pilon
```

```bash
# 步骤 1: BWA 比对
bwa mem -t 16 assembly.fasta \
  short_R1.fastq.gz short_R2.fastq.gz \
  > illumina_aligned.sam

# 步骤 2: SAMtools 排序索引
samtools view -bS illumina_aligned.sam | \
  samtools sort -@ 8 -o illumina_sorted.bam
samtools index illumina_sorted.bam

# 步骤 3: Pilon polishing
java -Xmx128G -jar pilon.jar \
  --genome assembly.fasta \
  --frags illumina_sorted.bam \
  --output pilon_polished \
  --changes \
  --vcf \
  --threads 16
```

**关键参数**:
| 参数 | 说明 |
|------|------|
| `--genome` | 输入组装 FASTA |
| `--frags` | Illumina 双端 BAM (插入片段 < 1000bp) |
| `--jumps` | 大插入片段文库 BAM |
| `--output` | 输出前缀 |
| `--changes` | 输出修改列表 |
| `--vcf` | 输出 VCF 格式的变异 |
| `--fix all` | 修复所有错误类型 |

### 4. Racon (快速备选)

**描述**: Racon 是超快的 polishing 工具，适合快速获取初步抛光结果。

```bash
racon -t 16 \
  long_reads.fastq.gz \
  overlaps.paf \
  assembly.fasta \
  > racon_polished.fasta
```

---

## Polishing 策略选择

| 数据场景 | 第一轮 | 第二轮 | 第三轮 |
|---------|--------|--------|--------|
| ONT 纯 | Flye (内置 1 round) | Medaka | Pilon (Illumina) |
| PacBio HiFi 纯 | hifiasm (无需) | gcpp | Pilon (可选) |
| ONT + Illumina | Flye | Medaka | Pilon |
| HiFi + Illumina | hifiasm | Pilon | - |
| ONT + HiFi | Flye | Medaka | Pilon |

### 多轮 Polishing 的必要性

每款工具修正的是不同类型的错误：
- **Medaka**: 修正 ONT 的随机单碱基错误（尤其在 homopolymer 区域）
- **gcpp**: 修正 PacBio HiFi 中残留的低质量位点
- **Pilon**: 用高准确度的短读长修正 SNV 和短 indel
- **Racon**: 快速修正长读长比对中的 mismatch

**建议**: 至少一轮 polishing 后检查 Merqury QV 和 BUSCO。如果 QV 已经 > 40，则不需要额外轮次。

---

## 植物基因组 Polishing 注意事项

1. **多倍体**: 多倍体会产生多个等位基因版本的 reads。Polishing 时需要保留杂合变异位点（不要过度平滑），设置合适的 `--mindepth` 和 `--minmq` 参数
2. **重复序列**: Pilon 在重复区域的变异调用可能有偏差（Illumina reads 比对到重复区域不可靠）
3. **叶绿体和线粒体**: 细胞器基因组通常是多拷贝，在 polishing 中会出现非常高的覆盖度
4. **着丝粒/端粒**: 这些区域的 reads 比对质量通常很低，polishing 效果有限

---

## Polishing 前后评估

```bash
# 使用 Merqury 评估 polishing 前后质量
merqury.sh genome.meryl assembly_raw.fasta raw_merqury
merqury.sh genome.meryl assembly_polished.fasta polished_merqury

# 比较 QV 值
grep "QV" raw_merqury.qv
grep "QV" polished_merqury.qv
```
