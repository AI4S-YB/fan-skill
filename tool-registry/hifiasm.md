# hifiasm -- 工具目录

## 概述

hifiasm 是 PacBio HiFi 长读长数据的专用组装器，专为高精度 HiFi reads（平均准确率 > 99.9%）设计。hifiasm 能在组装过程中同时解析单倍型（haplotype-resolved assembly），这对杂合度高的植物基因组特别有帮助。

## 推荐工具

### 1. hifiasm -- PacBio HiFi 组装金标准

**描述**: hifiasm 将 reads 的 overlap 关系编码为 phased string graph，在组装的同时分出两个单倍型。相较于传统的"组装后再分型"策略，hifiasm 的内置 phasing 更能保留杂合变异的真实性。

**安装**:
```bash
# 从源码编译
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make

# 或通过 conda
conda install -c bioconda hifiasm
```

**标准组装**:
```bash
# 基础 hifiasm 组装 (二倍体，无 Hi-C)
hifiasm -o assembly_prefix -t 32 \
  hifi_reads.fastq.gz

# 提取 primary assembly (bp.p_ctg)
awk '/^S/{print ">"$2;print $3}' \
  assembly_prefix.bp.p_ctg.gfa > assembly_primary.fasta
```

**各输出文件**:
| 文件 | 说明 | 推荐使用 |
|------|------|---------|
| `*.bp.p_ctg.gfa` | Primary contig graph | 是（primary assembly） |
| `*.bp.a_ctg.gfa` | Alternate contig graph | 备选（alternative assembly） |
| `*.bp.hap1.p_ctg.gfa` | Haplotype 1 resolved | phase-specific |
| `*.bp.hap2.p_ctg.gfa` | Haplotype 2 resolved | phase-specific |

**含 Hi-C 数据的组装**:
```bash
hifiasm -o assembly_prefix -t 32 \
  --h1 hic_R1.fastq.gz --h2 hic_R2.fastq.gz \
  hifi_reads.fastq.gz
```

Hi-C 数据可以帮助 hifiasm 区分同源染色体（haplotype phasing）并实现更长的连续性。

**含 HiFi + ONT 超长读长的组装**:
```bash
hifiasm -o assembly_prefix -t 32 \
  --ul ont_ultra_long.fastq.gz \
  hifi_reads.fastq.gz
```

### 关键参数

| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `-o` | 输出前缀 | 物种名缩写 |
| `-t` | 线程数 | 16-64 |
| `--h1/--h2` | Hi-C reads (R1/R2) | 可选 |
| `--ul` | 超长读长 (ONT > 100kb) | 可选，用于跨度 gap |
| `--n-hap` | 单倍型数 | 2 (二倍体), 4 (四倍体) |
| `-l` | 需保留的 reads 长度阈值 | 0 (保留所有) |
| `--hom-cov` | 同源覆盖倍数 | 自动 (通常不需设置) |
| `--purge-max` | Purge 最大深度 | 需时设置 |

### Purge Duplicates (去冗余)

对于高杂合度植物（杂合度 > 1%），hifiasm 可能会产生多套 contig。使用 purge_dups 去除冗余：

```bash
# 步骤 1: 使用 minimap2 比对
minimap2 -t 16 -x map-hifi \
  assembly_primary.fasta hifi_reads.fastq.gz \
  > alignment.paf

# 步骤 2: 计算覆盖直方图
calcuts coverage_histogram.txt < alignment.paf

# 步骤 3: 运行 purge_dups
purge_dups -T cutoffs -c coverage_histogram.txt \
  assembly_primary.fasta
```

---

## 植物 hifiasm 组装建议

### 基因组特征指导

| 植物基因组特征 | hifiasm 建议 |
|--------------|-------------|
| 低杂合度 (< 0.5%) | 标准参数，primary assembly 即可 |
| 中等杂合度 (0.5-1%) | 标准参数，检查 phased assemblies |
| 高杂合度 (> 1%) | 需 purge_dups，设置 `--purge-max` |
| 同源四倍体 | `--n-hap 4` |
| 异源多倍体 | 标准参数后亚基因组分离 |

### 测序深度建议

| 基因组大小 | 推荐 HiFi 覆盖度 |
|-----------|----------------|
| < 500 Mb | 20-30x |
| 500 Mb - 1 Gb | 20-30x |
| 1-3 Gb | 30-40x |
| > 3 Gb | 30-50x |

### 内存和运行时间估计

| 基因组大小 | 推荐 RAM | 预估运行时间 (40 线程) |
|-----------|---------|---------------------|
| 拟南芥 (135Mb) | 8 GB | 1-2 小时 |
| 水稻 (430Mb) | 16 GB | 4-8 小时 |
| 玉米 (2.4Gb) | 64 GB | 24-48 小时 |
| 小麦 (17Gb) | 256 GB+ | 5-10 天 |

---

## 常见问题

### Q: hifiasm "Killed" (OOM 错误)?
A: 内存不足。hifiasm 的 RAM 使用与基因组大小和杂合度成正比。对于 1Gb 以上的植物基因组，至少需要 64 GB RAM。可尝试使用 `-s` 参数降低内存使用（会牺牲部分组装质量）。

### Q: 产生大量小 contig (< 10kb) ?
A: 可能原因：(1) 测序深度不足，低于 20x；(2) reads 质量不好或不均一；(3) 物种杂合度过高。尝试提高测序覆盖度或使用 `--ul` 加入 ONT 超长读长。

### Q: Primary 与 Alternate 组装的取舍?
A: Primary assembly 包含更完整的基因组（但有部分 collapse），使用于大多数下游分析。Alternate assembly 保留了在 primary 中被 collapse 的区域。如果关注杂合变异，可以同时保留两套。
