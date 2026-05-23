# 亚硫酸盐测序数据预处理 -- 工具目录

## 概述

亚硫酸盐测序（WGBS/RRBS）数据的预处理包括原始数据质控、接头修剪、bisulfite 特异性的质量过滤。由于 bisulfite 转化改变了序列的碱基组成（未甲基化 C 转为 T），质控标准需要做特殊调整。

## 推荐流程

### 步骤 1: 原始数据质控

```bash
# FastQC 评估原始数据质量
fastqc -t 8 raw/*.fastq.gz -o qc/raw/
multiqc qc/raw/ -o qc/raw_report/
```

**检查要点（bisulfite 数据特殊考虑）**:
- **Per base sequence content**: bisulfite 转化的数据在 Read 1 的 5' 端会有 C 含量下降（因为非甲基化 C 被读成 T），这是正常现象而非质量问题
- **GC content**: WGBS 的 GC 分布可能显示为双峰（C->T 转化后），属于正常现象
- **Adapter content**: 仍需检查接头污染

### 步骤 2: 接头与质量修剪 (Trim Galore)

**Trim Galore** 是 bisulfite 数据修剪的专用工具，推荐使用：

```bash
# 双端 WGBS 修剪
trim_galore --paired \
  --clip_R1 10 --clip_R2 10 \
  --three_prime_clip_R1 10 --three_prime_clip_R2 10 \
  --fastqc \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  -o trimmed/

# RRBS 修剪（需要特殊处理 MspI 酶切位点）
trim_galore --paired \
  --rrbs \
  --fastqc \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  -o trimmed_rrbs/
```

**参数说明**:
| 参数 | 说明 |
|------|------|
| `--paired` | 双端模式 |
| `--clip_R1 N` | 切除 Read 1 5' 端 N bp（推荐 10bp） |
| `--three_prime_clip_R1 N` | 切除 Read 1 3' 端 N bp |
| `--rrbs` | RRBS 专用模式（去除 MspI 酶切位点引入的序列） |
| `--quality N` | 质量阈值（默认 20） |
| `--length N` | 最小 read 长度（默认 20） |
| `--fastqc` | 修剪后运行 FastQC |

### 步骤 3: 修剪后质控

```bash
# 再次运行 FastQC 验证修剪效果
fastqc -t 8 trimmed/*.fq.gz -o qc/trimmed/
multiqc qc/trimmed/ -o qc/trimmed_report/
```

### 步骤 4: 比对前统计

```bash
# 统计 reads 数
for f in trimmed/*.fq.gz; do
  count=$(zcat $f | wc -l)
  echo "$f: $((count / 4)) reads"
done
```

---

## WGBS vs RRBS 预处理差异

| 特性 | WGBS | RRBS |
|------|------|------|
| 覆盖范围 | 全基因组 | CpG 富集区域 |
| 文库大小 | 大（通常 >200M reads） | 较小（20-50M reads） |
| 对照样本 | 通常不需要 | 通常不需要 |
| 修剪工具 | Trim Galore (标准) | Trim Galore (--rrbs) |
| 比对方向 | non_directional（植物） | directional（取决于 protocol） |

---

## 植物 WGBS 特有考虑

### 亚硫酸盐转化效率评估

植物细胞壁较厚，DNA 提取过程中的杂质可能影响 bisulfite 转化效率。

**转化率检查** -- 叶绿体 DNA 作为内源对照：
植物叶绿体 DNA 几乎没有甲基化（或甲基化水平极低）。使用叶绿体基因组作为内源非甲基化对照来估算转化率：

```bash
# 比对到叶绿体参考基因组
bismark --bowtie2 -p 4 \
  --genome chloroplast_genome/ \
  sample_trimmed.fq.gz

# 叶绿体中残余的 C 甲基化应 < 1%
# 如果 > 1%，说明转化不完全
```

### 多倍体植物的覆盖度要求

植物多倍体（如小麦 6x）需要更深的测序深度：

| 倍性 | 物种示例 | 推荐深度 (WGBS) |
|------|---------|-----------------|
| 二倍体 | 拟南芥、番茄 | 15-30x |
| 二倍体 | 水稻、玉米（古四倍体但像二倍体） | 15-30x |
| 四倍体 | 棉花、油菜 | 30-50x |
| 六倍体 | 小麦 | 50-80x |

### 组织特定考虑

植物不同组织的甲基化谱差异显著：

- **叶片**: 最常见的组织，叶绿体 DNA 含量高（绿色叶片 ~30-50% reads 来自 cpDNA，需移除）
- **种子**: 胚与胚乳的甲基化模式完全不同
- **根**: 微生物污染风险（根际微生物 DNA）
- **花**: 生殖组织的甲基化与体细胞不同（涉及 imprinting）

---

## 质量指标

| 指标 | 优秀 | 可接受 | 需要关注 |
|------|------|--------|---------|
| Bismark 唯一比对率 | > 60% | 40-60% | < 40% |
| Bisulfite 转化率 | > 99% | 97-99% | < 97% |
| CG 甲基化 (叶片) | 20-50% | 10-60% | 异常极端 |
| Reads after trimming | > 95% | 85-95% | < 85% |
