# Bismark 比对与甲基化提取 -- 工具目录

## 概述

Bismark 是亚硫酸盐测序（WGBS/RRBS）数据比对和甲基化提取的金标准工具。它将 reads 同时比对到 bisulfite 转化后的参考基因组（C->T 和 G->A 两个版本），判定每个比对的最佳方向。

## 推荐工具

### 1. Bismark -- 金标准

**描述**: Bismark 将 reads 同时比对到参考基因组的 bisulfite 转换版本，通过多个比对的比较来确定正确的基因组位置和甲基化状态。

**安装**:
```bash
conda install -c bioconda bismark bowtie2 samtools
```

### 流程 1: WGBS 分析

#### 步骤 1: 基因组 Bisulfite 索引

```bash
# 构建 bisulfite 转换后的基因组索引
bismark_genome_preparation --path_to_aligner bowtie2 \
  --bowtie2 --verbose \
  genome_dir/
```

#### 步骤 2: Bismark 比对

```bash
# 双端 WGBS 比对
bismark --bowtie2 -p 8 \
  --non_directional \
  --genome genome_dir/ \
  -1 sample_R1_trimmed.fq.gz \
  -2 sample_R2_trimmed.fq.gz \
  --output_dir bismark_output/ \
  --temp_dir temp/

# 单端 WGBS 比对
bismark --bowtie2 -p 8 \
  --non_directional \
  --genome genome_dir/ \
  sample_trimmed.fq.gz \
  --output_dir bismark_output/
```

**关键参数**:

| 参数 | 说明 |
|------|------|
| `--bowtie2` | 使用 Bowtie2 作为比对器 |
| `-p 8` | 8 线程并行 |
| `--non_directional` | 非方向性文库（植物 WGBS 必须使用） |
| `--directional` | 方向性文库（如 TruSeq 甲基化建库） |
| `--pbat` | PBAT 文库模式 |
| `-I 0 -X 1000` | 插入片段长度范围 |
| `--score_min L,0,-0.6` | 比对得分阈值 |
| `--un` | 输出未比对 reads |
| `--ambiguous` | 输出多重比对 reads |

#### 步骤 3: 去重复

```bash
# 去 PCR 重复
deduplicate_bismark -p --bam \
  --output_dir bismark_output/ \
  bismark_output/sample_bismark_bt2_pe.bam
```

**参数**:
- `-p`: 双端模式
- `--bam`: 输出 BAM 格式
- `-s`: 单端模式

#### 步骤 4: 甲基化提取

```bash
# 提取甲基化信息
bismark_methylation_extractor \
  --paired-end --no_overlap \
  --comprehensive \
  --bedGraph --CX_context \
  --cytosine_report \
  --genome_folder genome_dir/ \
  --output_dir methylation_output/ \
  bismark_output/sample_deduplicated.bam
```

**关键参数**:
| 参数 | 说明 |
|------|------|
| `--paired-end` | 双端模式 |
| `--no_overlap` | 不重复计数双端重叠区 |
| `--comprehensive` | 输出所有甲基化信息 |
| `--bedGraph` | 输出 bedGraph 格式 |
| `--CX_context` | 输出所有上下文 (CG/CHG/CHH) |
| `--cytosine_report` | 生成 cytosine report |
| `--merge_non_CpG` | 合并非 CpG 甲基化 |

#### 步骤 5: Coverage 报告

```bash
# 生成 coverage 文件
coverage2cytosine --genome_folder genome_dir/ \
  -o sample_CX_report.txt \
  methylation_output/sample_deduplicated.bismark.cov.gz
```

### 流程 2: RRBS 分析

RRBS 使用 MspI 酶切富集 CpG 富集区域，比对时需要特殊处理酶切末端。

```bash
# RRBS 模式比对
bismark --bowtie2 -p 8 \
  --genome genome_dir/ \
  -1 sample_R1_trimmed.fq.gz \
  -2 sample_R2_trimmed.fq.gz \
  --rrbs \
  --output_dir bismark_rrbs_output/

# 后续步骤与 WGBS 相同，使用 bismark_methylation_extractor
```

---

## 植物甲基化 Bismark 特殊参数

### --non_directional (必须)

植物甲基化是非对称的（non-directional），意味着甲基化的 C 可以出现在 read1 或 read2 的任何位置。**植物 WGBS 几乎总是需要使用 `--non_directional` 参数**。

```bash
# 非方向性比对 = 4 个比对线程（每个 read 进行 OT 和 OB 比对）
# OT = Original Top strand, OB = Original Bottom strand
# CTOT = Complementary to OT, CTOB = Complementary to OB
bismark --non_directional ...
```

**注意**: 使用 `--non_directional` 会导致比对数为 4 倍，比对速度下降约 4 倍。

### 植物基因组大小参数

植物基因组大小变化极大（拟南芥 135Mb 到小麦 17Gb），这直接影响 Bismark 的 RAM 使用。

| 物种 | 基因组 | 推荐 RAM |
|------|--------|---------|
| 拟南芥 | 135 Mb | 4 GB |
| 水稻 | 430 Mb | 8 GB |
| 番茄 | 900 Mb | 16 GB |
| 玉米 | 2.4 Gb | 32 GB |
| 小麦 | 17 Gb | 128 GB+ |

---

## 甲基化上下文分离

植物 DNA 甲基化存在于三种序列上下文中。Bismark 可以输出这三种上下文的独立覆盖信息：

```bash
# 从 CX_report 分离各上下文
grep "CG_" sample_CX_report.txt > sample_CG.txt
grep "CHG_" sample_CX_report.txt > sample_CHG.txt
grep "CHH_" sample_CX_report.txt > sample_CHH.txt
```

---

## 比对报告解读

Bismark 会输出详细的比对报告。关键指标：

```
=== Bismark report ===
Total reads:                     50,000,000
Uniquely aligned reads:          22,500,000  (45.0%)
Reads with no alignment:         25,000,000  (50.0%)
Reads with multiple alignments:   2,500,000   (5.0%)

Cytosine methylation:
  Total Cs analysed:             500,000,000
  Methylated Cs in CpG:           80,000,000  (80.0%)
  Methylated Cs in CHG:           25,000,000  (25.0%)
  Methylated Cs in CHH:            5,000,000   (5.0%)
```

**质量基准**:
- **唯一比对率**: WGBS 40-60% 正常，RRBS 50-70% 预期
- **CG 甲基化**: 植物中 20-60%，因物种和组织而异
- **CHG 甲基化**: 植物中 10-30%
- **CHH 甲基化**: 植物中 1-10%
