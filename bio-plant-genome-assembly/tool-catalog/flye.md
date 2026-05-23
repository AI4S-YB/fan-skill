# Flye -- 工具目录

## 概述

Flye 是长读长（ONT 和 PacBio）基因组组装器，专为处理高错误率的长读长数据（ONT 错误率 5-15%）设计。Flye 使用 A-Bruijn graph 而非 overlap graph，能更好地解析重复序列结构。

## 推荐工具

### 1. Flye -- ONT 组装金标准

**描述**: Flye 通过迭代错误校正和 A-Bruijn 图谱构建来处理错误率较高的长读长数据。它首先使用长 k-mer 构建初步组装图，然后通过 repeat graph 解析来纠正重复序列诱导的错误。

**安装**:
```bash
conda install -c bioconda flye
```

**ONT 数据组装**:
```bash
# ONT 原始 reads (推荐)
flye --nano-raw ont_reads.fastq.gz \
  --genome-size 500m \
  --out-dir flye_output/ \
  --threads 32

# ONT 校正后的 reads
flye --nano-corr ont_corrected.fasta.gz \
  --genome-size 500m \
  --out-dir flye_output/ \
  --threads 32

# ONT HQ (高质量，如 Q20+)
flye --nano-hq ont_hq_reads.fastq.gz \
  --genome-size 500m \
  --out-dir flye_output/ \
  --threads 32
```

**PacBio 数据组装**:
```bash
# PacBio HiFi
flye --pacbio-hifi hifi_reads.fastq.gz \
  --genome-size 500m \
  --out-dir flye_output/ \
  --threads 32

# PacBio CLR (连续长读长)
flye --pacbio-raw pacbio_clr.fastq.gz \
  --genome-size 500m \
  --out-dir flye_output/ \
  --threads 32
```

### 关键参数

| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `--genome-size` | 估计基因组大小 | 带后缀 m/g（如 500m、2.4g） |
| `--threads` | CPU 线程数 | 16-64 |
| `--asm-coverage` | 目标组装覆盖度 | 自动（默认） |
| `--iterations` | Polish 迭代次数 | 1-2 |
| `--meta` | 宏基因组模式 | 不使用（植物基因组） |
| `--plasmid` | 质粒恢复模式 | 不使用 |
| `--scaffold` | 启用 scaffolding | 可选 |

**基因组大小参考**:
```bash
# 示例
flye --nano-raw ont_data.fastq.gz \
  --genome-size 1.19e8  # 也可以用 119m
```

### 输出文件

| 文件 | 说明 |
|------|------|
| `assembly.fasta` | 最终组装结果（推荐使用） |
| `assembly_graph.gfa` | 组装图（可用 Bandage 可视化） |
| `assembly_info.txt` | Contig 统计信息（长度、覆盖度、circular） |
| `params.json` | 运行参数记录 |

### 组装统计

```bash
# 查看组装统计
grep "Total\|N50\|Longest" flye.log

# 详细 contig 信息
head -20 assembly_info.txt
```

**assembly_info.txt 列说明**:
| 列 | 说明 |
|-----|------|
| #seq_name | Contig 名 |
| length | 长度 (bp) |
| cov. | k-mer 覆盖度 |
| circ. | 是否为环状 (Y/N) |
| repeat | 重复状态 (Y/N) |
| mult. | 多重度 |
| graph_path | 图路径 |

### 可选: 使用短读长 Polishing

Flye 内置了一轮 consensus polishing，但可以使用额外的短读长或 HiFi reads 进行后续 polishing：

```bash
# Polishing with Medaka (ONT)
medaka_consensus -i ont_reads.fastq.gz \
  -d flye_output/assembly.fasta \
  -o medaka_polished/ \
  -t 16

# Polishing with Pilon (Illumina)
bwa mem assembly.fasta short_R1.fq.gz short_R2.fq.gz \
  | samtools sort -o illumina_aligned.bam
java -jar pilon.jar --genome assembly.fasta \
  --frags illumina_aligned.bam \
  --output pilon_polished/
```

---

## 植物 Flye 组装建议

### 数据选择

| 数据类型 | Flye 模式 | 适用场景 |
|---------|----------|---------|
| ONT Guppy5+ HAC | `--nano-hq` | 最新 ONT 碱基识别 |
| ONT Guppy4- | `--nano-raw` | 老旧 ONT 数据 |
| ONT Q20+ 试剂 | `--nano-hq` | 新试剂盒，高准确率 |
| PacBio HiFi | `--pacbio-hifi` | 备选（优先使用 hifiasm） |
| PacBio CLR | `--pacbio-raw` | 老 PacBio 数据 |

### 估计基因组大小

使用 GenomeScope 或其他 k-mer 方法准确估计基因组大小。Flye 的 `--genome-size` 参数只需要粗略估计即可（不影响最终结果，只影响内存使用）。

**快速 k-mer 估计**:
```bash
# 使用 Jellyfish 快速估计
jellyfish count -C -m 21 -s 1G -t 16 \
  -o genome.jf short_reads.fq.gz
jellyfish histo genome.jf > genome.histo
# 上传到 GenomeScope 获得基因组大小
```

### 内存和运行时间

| 基因组大小 | 推荐 RAM | 运行时间 (32 threads, ONT) |
|-----------|---------|--------------------------|
| < 100 Mb | 8 GB | 2-4 小时 |
| 100-500 Mb | 16 GB | 6-12 小时 |
| 500 Mb - 1 Gb | 32 GB | 12-24 小时 |
| 1-3 Gb | 64 GB | 24-72 小时 |
| > 3 Gb | 128 GB+ | 3-10 天 |

---

## 常见问题

### Q: Flye 组装后出现大量"loop"和"repeat"标记？
A: 植物基因组的重复序列（TE、着丝粒、rDNA）在重复图中产生复杂结构。`circ.=N` 且 `repeat=Y` 的 contig 可能代表未完全解析的重复序列。可以保留，但标注为"低置信度"。

### Q: 组装结果大小远小于预期基因组大小？
A: 原因包括：(1) 测序覆盖度不够覆盖重复序列（植物中 TE 极多）；(2) `--asm-coverage` 参数默认截断了高覆盖度区域；(3) 高度重复的着丝粒、端粒、rDNA 区域在长读长组装中仍然困难。

### Q: Flye 报错 "not enough reads" 或 "coverage too low"?
A: 提高测序深度。ONT 数据推荐 30-60x 原始覆盖度。如果是低覆盖度 (<20x)，降低 `--asm-coverage` 或使用 HiFi 补充。

### Q: Flye vs hifiasm vs Canu?
A: 选择参考：
- **ONT 数据**: Flye（最优）
- **PacBio HiFi 数据**: hifiasm（最优）
- **PacBio CLR 数据**: Flye 或 Canu
- **低覆盖度 + 高错误率**: Flye（最稳）
- **高杂合度二倍体**: hifiasm（phasing 更好）
