# 转座子分析 -- 分析笔记本

## 分析概览

本笔记本提供植物转座子(TE)分析和注释的完整流程，涵盖TE鉴定、分类、注释、插入时间分析和分布分析。

---

## 1. TE从头鉴定

### 1.1 EDTA流程 (推荐用于植物)

```bash
# EDTA完整流程
EDTA.pl --genome genome.fasta \
  --species plant \
  --sensitive 1 \
  --overwrite 1 \
  --anno 1 \
  --cds gene_cds.fasta \
  --threads 16

# 输出文件：
# genome.fasta.EDTA.TElib.fa - TE库
# genome.fasta.EDTA.TEanno.gff3 - TE注释
# genome.fasta.mod.EDTA.RM.out - RepeatMasker结果
```

### 1.2 RepeatModeler

```bash
# 构建数据库
BuildDatabase -name genome_db genome.fasta

# 运行RepeatModeler
RepeatModeler -database genome_db \
  -pa 16 \
  -LTRStruct \
  -gff

# 输出：
# genome_db-families.fa - TE库
# genome_db-families.stk - 多序列比对
```

### 1.3 特定类型TE鉴定

```bash
# LTR逆转录转座子
LTRharvest -index genome_index \
  -seq genome.fasta \
  -out ltr_output.gff \
  -outinner ltr_inner.fa

# LTR_retriever过滤
LTR_retriever -genome genome.fasta \
  -inharvest ltr_output.gff \
  -outgenome genome_ltr.fa

# MITE分析
MITE-Hunter.pl -genome genome.fasta -p 16

# Helitron分析
helitron_scanner.pl -genome genome.fasta -o helitron_output
```

---

## 2. TE分类注释

### 2.1 DeepTE深度学习分类

```bash
# DeepTE分类
DeepTE.py -i te_sequences.fa \
  -o deepte_output \
  -sp P

# 参数说明：
# -sp P: 植物
# -sp A: 动物
# -sp F: 真菌
```

### 2.2 RepeatClassifier

```bash
# 使用RepeatClassifier分类
RepeatClassifier -consensi te_consensi.fa \
  -rfam_abundant \
  -db reposted RepeatMaskerLib.h5
```

### 2.3 TE分类系统

| Class | Order | Superfamily | 示例 |
|-------|-------|-------------|------|
| I (逆转录) | LTR | Copia, Gypsy | Tnt1, Karma |
| I | non-LTR | LINE, SINE | LINE-1 |
| II (DNA) | TIR | hAT, Mutator, CACTA | Ac/Ds |
| II | non-TIR | Helitron, Maverick | Helitron |

---

## 3. 基因组TE注释

### 3.1 RepeatMasker

```bash
# 使用自定义TE库注释
RepeatMasker -pa 16 \
  -lib custom_te_lib.fa \
  -gff \
  -dir rm_output \
  genome.fasta

# 输出文件：
# genome.fasta.out - 详细注释
# genome.fasta.gff - GFF格式
# genome.fasta.masked - masked基因组
```

### 3.2 注释结果统计

```bash
# 统计TE覆盖度
python calculate_te_coverage.py genome.fasta.out

# 按类型统计
awk '{sum[$11]+=$7-$6} END{for(t in sum) print t, sum[t]}' genome.fasta.out | \
  sort -k2 -nr > te_stats.txt
```

---

## 4. LTR插入时间分析

### 4.1 提取完整LTR

```bash
# 从EDTA或LTR_retriever结果提取完整LTR
grep "complete" ltr_retriever.gff | \
  cut -f1,4,5,9 > complete_ltr.txt
```

### 4.2 计算插入时间

```bash
# 提取LTR序列
bedtools getfasta -fi genome.fasta \
  -bed complete_ltr.txt \
  -fo complete_ltr.fa

# 计算LTR对的序列分化度
python calculate_ltr_divergence.py complete_ltr.fa

# 转换为插入时间
# Time = K / (2 * substitution_rate)
# 植物常用替代率: 1.3e-8
python estimate_insertion_time.py \
  --divergence ltr_divergence.txt \
  --rate 1.3e-8 \
  --output insertion_times.txt
```

### 4.3 插入时间可视化

```R
# R绘图
library(ggplot2)
insertion <- read.table("insertion_times.txt", header=TRUE)
ggplot(insertion, aes(x=time_mya)) +
  geom_histogram(binwidth=0.1, fill="steelblue") +
  theme_minimal() +
  labs(x="Insertion Time (Mya)", y="LTR Count")
```

---

## 5. TE分布分析

### 5.1 染色体分布

```bash
# 滑动窗口统计TE密度
bedtools makewindows -g chromosome.sizes -w 100000 > windows.bed
bedtools coverage -a windows.bed -b te_annotation.gff > te_density.txt

# 可视化
Rscript plot_te_distribution.R te_density.txt
```

### 5.2 基因邻近分析

```bash
# 分析基因上下游的TE分布
# 提取基因位置
awk '$3=="gene"' genes.gff > gene_positions.bed

# 计算基因到最近TE的距离
bedtools closest -a gene_positions.bed -b te_annotation.bed -d > gene_te_distance.txt

# 统计分析
Rscript analyze_gene_te_distance.R gene_te_distance.txt
```

### 5.3 着丝粒区域TE

```bash
# 提取着丝粒区域TE
bedtools intersect -a te_annotation.bed -b centromere.bed > centromere_te.bed

# 统计着丝粒TE组成
awk '{sum[$4]++} END{for(t in sum) print t, sum[t]}' centromere_te.bed
```

---

## 6. TE-基因关联分析

### 6.1 TE来源的调控元件

```bash
# 提取启动子区域的TE
# 启动子定义：基因上游2kb
bedtools flank -g chromosome.sizes -b 2000 -i gene_positions.bed > promoters.bed
bedtools intersect -a promoters.bed -b te_annotation.bed > te_in_promoters.bed

# 统计
wc -l te_in_promoters.bed
```

### 6.2 TE表达分析

```bash
# 如果有RNA-seq数据
# 使用TE专用定量工具
salmon index -t te_library.fa -i te_index
salmon quant -i te_index -l A -r rnaseq_R1.fastq.gz -o te_quant

# 或使用TEtranscripts
TEtranscripts --BAM rnaseq.bam --TE te_annotation.gtf --project te_expression
```

### 6.3 TE甲基化分析

```bash
# 如果有BS-seq数据
# 提取TE区域甲基化
bedtools intersect -a methylation.bed -b te_annotation.bed -wa -wb > te_methylation.txt

# 计算平均甲基化水平
awk '{sum+=$4; count++} END{print sum/count}' te_methylation.txt
```

---

## 7. TE景观图

### 7.1 Kimura距离分析

```bash
# RepeatMasker生成Kimura距离
RepeatMasker -pa 16 -lib te_lib.fa -a genome.fasta
calcDivergenceFromAlign.pl -a genome.fasta.align > kimura_divergence.txt

# 生成景观图
rmOut2plot.R genome.fasta.out

# 或使用自定义脚本
Rscript plot_te_landscape.R kimura_divergence.txt te_classification.txt
```

---

## 8. 植物特异性分析

### 8.1 高重复基因组处理

```bash
# 对于高重复基因组（>50% TE）
# 增加内存和运行时间

# EDTA大基因组模式
EDTA.pl --genome large_genome.fa \
  --sensitive 0 \
  --repeatmaker 1 \
  --threads 32 \
  --maxln 20000
```

### 8.2 近期TE爆发分析

```bash
# 识别近期活跃的TE
# 基于低分化和完整结构

# 筛选完整LTR
grep "intact" ltr_retriever.gff | \
  awk '$9 < 0.01' > recent_ltr.gff  # 低分化
```

---

## 常见问题

### Q: TE注释比例过低？
A: 检查TE库质量。使用EDTA从头鉴定，或添加同源物种的TE库。

### Q: LTR插入时间分布异常？
A: 检查替代率是否合适。植物不同物种替代率可能有差异。

### Q: TE分类结果混乱？
A: 使用DeepTE等工具进行深度学习分类，提高准确率。

### Q: 重复序列比例异常高？
A: 检查是否有低复杂序列污染。使用dustmasker过滤简单重复。

---

## 参考

- EDTA: https://github.com/oushujun/EDTA
- RepeatModeler: https://www.repeatmasker.org/RepeatModeler/
- RepeatMasker: https://www.repeatmasker.org/
- DeepTE: https://github.com/LiLabAtVT/DeepTE
- LTR_retriever: https://github.com/oushujun/LTR_retriever
