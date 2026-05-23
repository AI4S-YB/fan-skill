# 植物基因组组装 -- 分析笔记本

## 分析概览

本笔记本提供植物基因组组装的完整流程，涵盖三代测序数据（PacBio HiFi 和 ONT）的从头组装、混合组装、碱基校正（polishing）、Hi-C 染色体挂载和质量评估。

---

## 1. 测序数据预处理

### 1.1 数据质控

```bash
# PacBio HiFi 数据
# HiFi reads 已经是高质量 CCS reads，通常不需要额外修剪
samtools fastq pacbio_hifi.bam > hifi_reads.fastq

# ONT 数据质控
# NanoPlot 评估 reads 质量和长度分布
NanoPlot -t 8 --fastq ont_reads.fastq.gz \
  -o nanoplot_output/ \
  --loglength --N50

# 过滤短 reads（可选）
seqkit seq -m 10000 ont_reads.fastq.gz > ont_reads_filt.fastq
```

### 1.2 短读长数据质控（用于混合组装）

```bash
# Trimmomatic 修剪
trimmomatic PE -threads 8 \
  short_R1.fq.gz short_R2.fq.gz \
  short_R1_trimmed.fq.gz short_R1_unpaired.fq.gz \
  short_R2_trimmed.fq.gz short_R2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 1.3 基因组调查 (Genome Survey)

在组装前执行 genome survey 获得 k-mer 频谱和基因组特征估计：

```bash
# Jellyfish + GenomeScope 2.0
jellyfish count -C -m 21 -s 1G -t 16 \
  -o genome.jf short_reads.fq.gz

jellyfish histo -t 16 genome.jf > genome.histo

# 上传 genome.histo 到 GenomeScope 在线工具：
# http://qb.cshl.edu/genomescope/genomescope2.0/
```

**GenomeScope 产出信息**:
- 估计基因组大小 (Estimated genome size)
- 杂合度 (Heterozygosity)
- 重复序列比例 (Repeat content)
- 倍性推断 (Ploidy)

---

## 2. 基因组组装

### 2.1 hifiasm (PacBio HiFi 最优)

hifiasm 是 PacBio HiFi 数据组装的推荐工具，特别适合高杂合度植物基因组。hifiasm 内置了 phased assembly，可直接产生 haplotype-resolved 组装结果。

```bash
# 标准 hifiasm 组装
hifiasm -o assembly_prefix -t 32 \
  hifi_reads.fastq.gz

# 含 Hi-C 数据的 hifiasm (Hi-C 辅助 phasing)
hifiasm -o assembly_prefix -t 32 \
  --h1 hic_R1.fq.gz --h2 hic_R2.fq.gz \
  hifi_reads.fastq.gz

# 提取 primary assembly
awk '/^S/{print ">"$2;print $3}' \
  assembly_prefix.bp.p_ctg.gfa > assembly_primary.fasta
```

**hifiasm 参数说明**:
- `-t`: 线程数
- `-o`: 输出前缀
- `--h1 / --h2`: Hi-C 读段文件（可选）
- `--ul`: 加入超长 ONT 读段辅助（可选）
- `--n-hap`: phasing 分区数量，默认为 2（二倍体）

**输出文件**:
- `*.bp.p_ctg.gfa`: Primary contig graph（推荐使用）
- `*.bp.a_ctg.gfa`: Alternate contig graph
- `*.bp.hap*.p_ctg.gfa`: Haplotype-resolved contig graphs

### 2.2 Flye (ONT 最优)

Flye 专门针对 ONT 长读长的错误模式进行了优化。

```bash
# ONT 数据 Flye 组装
flye --nano-raw ont_reads.fastq.gz \
  --genome-size 500m \
  --out-dir flye_output/ \
  --threads 32

# PacBio HiFi 数据 Flye (备选)
flye --pacbio-hifi hifi_reads.fastq.gz \
  --genome-size 500m \
  --out-dir flye_output/ \
  --threads 32
```

**Flye 参数说明**:
- `--nano-raw`: ONT 原始 reads
- `--nano-corr`: 已校正 ONT reads
- `--pacbio-raw`: PacBio CLR reads
- `--pacbio-corr`: 已校正 PacBio reads
- `--pacbio-hifi`: PacBio HiFi reads
- `--genome-size`: 估计基因组大小（可用后缀 m/g，如 500m）
- `--asm-coverage`: 目标覆盖度（默认自动检测）
- `--scaffold`: 长距离 scaffolding（可选）

**输出文件**:
- `assembly.fasta`: 最终组装结果
- `assembly_graph.gfa`: 组装图
- `assembly_info.txt`: Contig 信息和覆盖度统计

### 2.3 混合组装 (长读长 + 短读长)

当同时拥有长读长和短读长数据时，混合组装可以结合长读长的连续性优势与短读长的准确性优势。

```bash
# MaSuRCA 混合组装器
masurca -g masurca_config.txt

# config.txt 内容示例：
# DATA
# PE= pe 300 50 short_R1.fq.gz short_R2.fq.gz
# PACBIO= hifi_reads.fastq.gz
# END
# PARAMETERS
# GRAPH_KMER_SIZE = auto
# USE_LINKING_MATES = 1
# LIMIT_JUMP_COVERAGE = 300
# CA_PARAMETERS = ovlMerSize=30
# KMER_COUNT_THRESHOLD = 1
# NUM_THREADS = 32
# JF_SIZE = 10G
# SOAP_ASSEMBLY = 0
# END
```

---

## 3. 组装后处理 (Polishing)

### 3.1 Medaka (ONT Polishing)

```bash
# ONT reads 的碱基校正
medaka_consensus -i ont_reads.fastq.gz \
  -d flye_assembly.fasta \
  -o medaka_output/ \
  -t 16 \
  -m r941_min_hac_g507
```

### 3.2 gcpp (PacBio HiFi Polishing)

```bash
# HiFi reads 的碱基校正
gcpp -j 16 \
  --algorithm arrow \
  -r hifi_reads.bam \
  -o polished_assembly.fasta \
  hifiasm_assembly.fasta
```

### 3.3 Pilon (短读长 Polishing)

```bash
# 使用 Illumina short reads 进行 polishing
# 步骤 1: 比对
bwa mem -t 16 assembly.fasta short_R1.fq.gz short_R2.fq.gz \
  > alignment.sam

samtools view -bS alignment.sam | samtools sort -o alignment_sorted.bam
samtools index alignment_sorted.bam

# 步骤 2: Pilon polishing
java -Xmx64G -jar pilon.jar \
  --genome assembly.fasta \
  --frags alignment_sorted.bam \
  --output pilon_output \
  --threads 16
```

---

## 4. Hi-C 染色体挂载 (Scaffolding)

### 4.1 Hi-C 数据预处理

```bash
# Hi-C reads 比对
bwa mem -t 16 -5SP assembly.fasta hic_R1.fq.gz hic_R2.fq.gz \
  > hic_alignment.sam

# 或使用 Juicer 完整流程
juicer.sh -t 16 \
  -g genome_name \
  -s MboI \
  -z assembly.fasta \
  -y restriction_sites.txt \
  -p chrom.sizes \
  -d juicer_work/
```

### 4.2 YAHS (Yet Another Hi-C Scaffolder)

```bash
# YAHS 染色体挂载
yahs assembly.fasta hic_alignment.bam \
  -o yahs_output \
  --no-contig-ec
```

### 4.3 SALSA2 (备选)

```bash
# SALSA2 Hi-C scaffolding
samtools index hic_alignment.bam

python run_pipeline.py \
  -a assembly.fasta \
  -l assembly.fasta.fai \
  -b hic_alignment.bam \
  -e GATC \
  -o salsa_output/
```

---

## 5. 组装质量评估

### 5.1 BUSCO 完整性评估

```bash
# 使用植物 lineage 数据库进行 BUSCO 评估
busco -i assembly.fasta \
  -l embryophyta_odb10 \
  -o busco_output \
  -m genome \
  -c 16
```

**常用植物 BUSCO lineage**:
| Lineage | 适用物种 |
|---------|---------|
| embryophyta_odb10 | 所有陆地植物 |
| eudicots_odb10 | 双子叶植物 |
| liliopsida_odb10 | 单子叶植物 |
| brassicales_odb10 | 十字花科 |
| poales_odb10 | 禾本科 |
| solanales_odb10 | 茄科 |
| fabales_odb10 | 豆科 |

**BUSCO 评估标准**:
- C (Complete): > 90% 为优秀
- C (Complete): > 80% 为良好
- Single (S): 越高越好（>70% 理想，但有重复是正常的）
- Duplicated (D): 对二倍体应低（<10%），多倍体中 D 升高是合理的
- Fragmented (F): < 10% 为理想
- Missing (M): < 10% 为理想

### 5.2 Merqury k-mer 评估

```bash
# 首先生成 k-mer 数据库
meryl k=21 count output genome.meryl short_reads.fq.gz

# Merqury 评估
merqury.sh genome.meryl assembly.fasta output_prefix
```

**Merqury 产出指标**:
- **QV (Consensus Quality Value)**: 组装碱基准确度，> 30 为较好，> 40 为优秀
- **k-mer 完整性 (Completeness)**: 表示组装覆盖了多少 reads 中的 k-mer，> 90% 为良好
- **k-mer 谱图 (Spectra-cn)**: 可视化倍性和重复

### 5.3 基本统计

```bash
# 组装统计
assembly-stats assembly.fasta > assembly_stats.txt

# N50, L50, 总长度, 最大 contig 等
cat assembly_stats.txt
```

**组装质量金标准**:
| 指标 | 优秀 | 良好 | 勉强可接受 |
|------|------|------|-----------|
| Contig N50 | > 10 Mb | > 1 Mb | > 100 kb |
| BUSCO Complete | > 95% | > 90% | > 80% |
| Merqury QV | > 40 | > 30 | > 20 |
| k-mer Completeness | > 95% | > 90% | > 80% |
| Scaffold N50 (若含 Hi-C) | > 50 Mb | > 20 Mb | > 10 Mb |

---

## 6. 可视化

### 6.1 组装图可视化

```bash
# Bandage 可视化组装图
Bandage image assembly_graph.gfa assembly_graph.png \
  --height 3000 --width 4000
```

### 6.2 Hi-C Contact Map

```bash
# HiCExplorer 生成 contact map
hicBuildMatrix -s hic_alignment.bam \
  --binSize 50000 \
  -o hic_matrix.h5

hicPlotMatrix -m hic_matrix.h5 \
  -o hic_contact_map.png \
  --log1p --dpi 300
```

---

## 常见问题

### Q: hifiasm 组装结果严重碎片化怎么办？
A: 检查测序深度。HiFi reads 推荐覆盖度 30-40x。如果覆盖度过低（< 20x），组装连续性会明显下降。检查物种杂合度——高杂合度植物基因组中 hifiasm 可能产生更多的 haplotype separation。

### Q: Flye 组装中出现"loops"（组装图中的环）？
A: 植物基因组中大量重复序列会导致组装图出现复杂结构。尝试调整 `--asm-coverage` 为 30-40x（而非默认的自动检测），以及在 scaffolding 后使用 purge_dups 去除冗余。

### Q: 组装大小远大于或远小于预期？
A: 大于预期：可能存在 phasing 导致的 haplotype duplication，使用 `purge_dups` 清理。或杂合度过高导致的两套单倍型被组装为独立 contig。小于预期：检查测序覆盖度是否足够覆盖高度重复区域。植物基因组中着丝粒和 rDNA 区域极难组装。

### Q: BUSCO Duplicated 比例过高？
A: 对于二倍体物种，如果 DUPLICATED > 20%，可能存在以下问题：
1. 杂合度极高导致 redundant haplotypes（使用 purge_dups）
2. 污染了不同基因型的 reads
3. 真有全基因组重复（WGD）历史——检查物种系统发育背景

### Q: Hi-C 挂载后出现大片段 misassembly？
A: 检查 Hi-C contact map。对角线外的强信号表示 misassembly。使用 Juicebox 手动校正。或使用 Manual curation 工具（如 Pretext, JBrowse2）。

---

## 参考

- hifiasm: https://github.com/chhylp123/hifiasm
- Flye: https://github.com/fenderglass/Flye
- Medaka: https://github.com/nanoporetech/medaka
- BUSCO: https://busco.ezlab.org/
- Merqury: https://github.com/marbl/merqury
- YAHS: https://github.com/c-zhou/yahs
- GenomeScope 2.0: http://qb.cshl.edu/genomescope/genomescope2.0/
- purge_dups: https://github.com/dfguan/purge_dups
