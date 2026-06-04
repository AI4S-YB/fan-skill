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

---

## 植物基因组组装的特有考量

### 非模式植物组装策略

大多数植物研究者面对的是非模式物种——没有近缘参考基因组，也没有现成的组装方案：

**基因组调查（Genome Survey）先行**：
- 在启动大规模测序前，先用 Illumina short reads（30-50x coverage）做 GenomeScope 分析
- 获取以下关键参数：估计基因组大小、杂合度、重复序列比例、倍性
- 基因组大小直接影响测序量和计算资源需求：< 1Gb 的基因组在单机上可完成，> 5Gb 需要集群
- 高杂合度（> 1%）的物种需要更高的长读长覆盖度（50-60x HiFi 而非 30-40x）来区分 haplotypes

**测序策略选择考虑因素**：

| 基因组特征 | 推荐策略 | 原因 |
|-----------|---------|------|
| 小基因组 (<500Mb) + 低杂合 | HiFi 30x | 单 SMRT cell 8M 即可完成 |
| 小基因组 + 高杂合 | HiFi 40-50x + Hi-C | 高杂合需要更多 reads 分型 |
| 中等基因组 (500Mb-2Gb) + 二倍体 | HiFi 30-40x + Hi-C | 标准植物基因组方案 |
| 大基因组 (2-5Gb) + 高重复 | HiFi 30x + ONT 20x 混合 | 互补长读长技术克服重复区域 |
| 超大基因组 (>5Gb) | ONT 50-60x | 成本可控（HiFi 在此规模成本过高） |
| 多倍体 (≥4x) | HiFi 50-60x + Hi-C + 亲本 Illumina | 亲本 reads 辅助亚基因组分型 |

**低预算替代方案**：
- 如果预算仅够 Illumina short reads：组装质量将显著受限（N50 在 Kb 级别而非 Mb）。仅建议用于基因空间的快速调查（gene space survey），不能发表为参考基因组。
- ONT MinION 是成本最低的长读长方案：FLOW-MIN114 (R10.4.1) 可产出 15-30Gb/flow cell，对小基因组（<500Mb）可以做到 50x 覆盖
- 合租测序平台（如多家实验室共享一个 PacBio Revio run）可大幅降低单物种成本

### 育种群体与品系选择的组装考量

植物基因组组装不同于动物——常使用育种品系或特定基因型：

**自交系 vs 杂交种**：
- 自交系（如玉米 B73、水稻 Nipponbare）：杂合度低，组装相对容易。是最理想的材料。
- 杂交种（F1）：杂合度高，hifiasm 会产生两个 haplotypes 的组装。如果目标是获得参考组装，选亲本自交系而非杂交种。
- 高度杂合的无性繁殖作物（马铃薯、甘薯、甘蔗）：无法获得纯合个体，只能接受 haplotype-resolved 或 collapsed assembly。

**DH 群体（双单倍体）**：
- 对于可产生 DH 的作物（玉米、油菜、大麦），DH 系是理想材料（纯合、无杂合度）
- 组装质量通常最高（N50 最高，BUSCO Duplicated 最低）

**异交物种（如黑麦草、苜蓿、向日葵野生种）**：
- 由于自交不亲和，每个个体都是高度杂合的
- 建议从群体中选择一个"代表个体"进行单个体深度测序
- 或者采用 trio-binning 策略（测序亲本 + F1 后代），利用亲本 k-mer 将 F1 的 reads 分箱后再组装

### 植物特有基因组区域的组装困难

某些植物基因组区域即使使用 HiFi + Hi-C 也难以完整解析：

**rDNA 区域**：
- 植物中 45S rDNA（18S-5.8S-25S）和 5S rDNA 串联重复数千拷贝
- 即使 HiFi reads（15-25Kb）也无法跨越完整的 rDNA array（通常 5-10Mb）
- 策略：接受 rDNA 区域的 gap；或用 ONT 超长 reads（>100Kb）尝试跨越
- 在 BUSCO 评估中排除 rDNA 相关的 BUSCO 基因

**着丝粒区域**：
- 植物着丝粒由卫星重复序列（CentC 等）和 CRM 逆转录转座子组成，长度可达数 Mb
- 完全解析着丝粒通常需要 ONT ultra-long reads (>100Kb) 或专门的长距离连锁方法
- 对于大多数非模式植物组装，着丝粒 gap 是可以接受的

**端粒区域**：
- 植物端粒重复序列为 TTTAGGG（拟南芥型）或 TTTTAGGG（某些单子叶植物型）
- hifiasm 和 Flye 通常在端粒处终止
- 可以通过识别 contig 末端的端粒重复来判断组装是否达到染色体末端

### 细胞器基因组组装

植物含有三个基因组：

- **叶绿体基因组**：通常 120-160Kb，环状。在 HiFi reads 中覆盖度极高（50-200x 于核基因组）。在组装核基因组前先组装叶绿体基因组（使用 OATK 或 GetOrganelle），因为 cpDNA reads 会干扰核基因组组装 graph。
- **线粒体基因组**：植物线粒体基因组很大（200Kb-2.7Mb）且结构复杂（多分支、重组），组装难度最高。通常需要专门的线粒体组装工具或混合组装策略。
- **核基因组**：主组装目标。

**工作流程建议**：先分离组装叶绿体和线粒体基因组，过滤掉细胞器 reads 后再进行核基因组组装。

### 作物特异性组装资源

| 资源 | 适用作物 | 内容 |
|------|---------|------|
| Plant Genome Portal (Phytozome) | 所有植物 | 300+ 已注释植物基因组下载 |
| Gramene | 禾本科 | 基因组、变异、表达数据库 |
| SoyBase | 大豆 | 大豆基因组和遗传学资源 |
| MaizeGDB | 玉米 | 玉米 B73 和泛基因组 |
| WheatIS | 小麦 | 小麦基因组资源整合 |
| Brassica Database (BRAD) | 芸苔属 | 油菜、白菜、甘蓝基因组 |
| Sol Genomics Network | 茄科 | 番茄、马铃薯、辣椒、茄子 |

---

## 多平台组装评估矩阵

对于高重要性基因组项目，建议运行多个组装工具并比较结果：

### 组装工具对比评估

```bash
# 运行多个组装工具
# 1. hifiasm (HiFi首选)
hifiasm -o hifiasm_asm -t 32 hifi_reads.fastq.gz

# 2. Flye (ONT首选，也可用于HiFi)
flye --pacbio-hifi hifi_reads.fastq.gz \
  --genome-size 500m \
  --out-dir flye_output \
  --threads 32

# 3. IPA (PacBio官方流程)
ipa local --fasta hifi_reads.fastq.gz --threads 32

# 4. NEAT (噪声感知组装)
neat assemble -i hifi_reads.fastq.gz -o neat_output
```

### 评估指标矩阵

| 指标 | hifiasm | Flye | IPA | 说明 |
|------|---------|------|-----|------|
| Contig N50 | 记录 | 记录 | 记录 | 越高越好 |
| BUSCO Complete | 记录 | 记录 | 记录 | >90%为良好 |
| Merqury QV | 记录 | 记录 | 记录 | >40为优秀 |
| 组装总长度 | 记录 | 记录 | 记录 | 与估计大小比较 |
| LTR回收率 | 记录 | 记录 | 记录 | 植物重要指标 |

### 结果整合策略

```bash
# 如果一个组装明显优于其他，直接选用
# 如果各有所长，考虑：

# 1. GFAtools合并
gfatools merge asm1.gfa asm2.gfa -o merged.gfa

# 2. 选择性拼接
# 对特定区域选择最佳组装片段

# 3. 差异报告
# 详细报告各工具优劣，由用户决定
```

---

## T2T (Telomere-to-Telomere) 组装流程

追求完整染色体级别的端粒到端粒组装：

### 1. Gap填补

```bash
# 使用ONT超长reads填补gap
# 首先识别gap位置
grep "N" assembly.fasta | head -20

# 使用TGS-GapCloser或LR_Gapcloser
LR_Gapcloser.sh -i assembly.fasta \
  -l ont_ultra_long.fastq.gz \
  -o gap_closed.fasta \
  -t 16

# 或使用TGS-GapCloser
python TGS-GapCloser.py \
  --scaff assembly.fasta \
  --reads ont_ultra_long.fastq.gz \
  --output gap_closed \
  --threads 16
```

### 2. 端粒验证

```bash
# 检测端粒重复序列
# 植物端粒序列：TTTAGGG (拟南芥型)
# 某些单子叶植物：TTTTAGGG

# 搜索contig末端的端粒重复
# 正向端粒
grep -E "^TTTAGGG{10,}" assembly.fasta

# 反向端粒
grep -E "CCCTAAA{10,}$" assembly.fasta

# 使用FindTelomere工具
find_telomeres.py -i assembly.fasta -r TTTAGGG -m 10 -o telomere_report.txt
```

### 3. 着丝粒区域分析

```bash
# 识别着丝粒卫星重复
# 使用TAREAN (Terminal-Associated REpeat ANalyzer)
tarean -i Illumina_reads.fq -o tarean_output

# 检测着丝粒特异重复
# 玉米：CentC, CRM1/CRM2
# 小麦：CCS1, 180bp重复
# 水稻：CentO

# 使用RepeatMasker注释着丝粒区域
RepeatMasker -species "Oryza sativa" \
  -pa 8 \
  assembly.fasta
```

### 4. T2T完整性评估

```bash
# 端粒数量统计（预期=2×染色体数）
expected_telomeres=$((2 * chromosome_count))
actual_telomeres=$(grep -c "telomere_signal" telomere_report.txt)
echo "端粒完整性: $actual_telomeres / $expected_telomeres"

# 着丝粒检测
# 预期每条染色体有1个着丝粒区域

# Gap统计
gap_count=$(grep -c "^>" assembly.fasta | xargs -I {} sh -c 'grep -o "N\+" assembly.fasta | wc -l')
echo "剩余gap数量: $gap_count"
```

---

## 多倍体组装策略

### 同源多倍体 vs 异源多倍体

| 类型 | 基因组结构 | 组装难度 | 推荐策略 |
|------|-----------|---------|---------|
| 同源多倍体 | AAAA (相同或相似) | 极高 | hifiasm --n-hap, 预期多单倍型 |
| 异源多倍体 | AABB (不同来源) | 中等 | 按亚基因组分开处理 |
| 节段异源多倍体 | AA BB (部分同源) | 高 | 结合两种策略 |

### 同源四倍体组装示例

```bash
# 同源四倍体（如马铃薯）
# hifiasm会尝试分出4个单倍型
hifiasm -o tetraploid_asm \
  --n-hap 4 \
  -t 32 \
  hifi_reads.fastq.gz

# 注意：结果可能非常破碎
# 替代策略：接受collapsed assembly

# 使用purge_dups控制冗余
purge_dups -T 0 -e 0.5 assembly.fasta
```

### 异源四倍体组装示例

```bash
# 异源四倍体（如棉花AADD、油菜AACC）
# 可以按亚基因组分开处理

# 方法1：直接组装后分离
hifiasm -o allo_asm -t 32 hifi_reads.fastq.gz
# 后续根据共线性和BUSCO分布分离亚基因组

# 方法2：Trio-binning（需要亲本数据）
# 测序两个亲本（如A基因组和D基因组来源）
# 根据亲本k-mer分离reads
trio_binning.py -p parentA.fq -m parentD.fq -c child.fq

# 分别组装两个bin
hifiasm -o hapA_asm -t 32 binA.fastq.gz
hifiasm -o hapD_asm -t 32 binD.fastq.gz
```

### 多倍体组装质量评估

```bash
# BUSCO评估多倍体
busco -i assembly.fasta \
  -l embryophyta_odb10 \
  -m genome \
  -c 16

# 注意：多倍体BUSCO Duplicated会很高，这是正常的
# 四倍体预期Duplicated: 40-60%

# 使用单拷贝核基因评估
# 选择植物单拷贝基因集，检查每个基因的拷贝数
```
| CottonGen | 棉花 | 棉花基因组和育种资源 |
| GrainGenes | 小谷物 | 大麦、燕麦、黑麦、小麦 |
