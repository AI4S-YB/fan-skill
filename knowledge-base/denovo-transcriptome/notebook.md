# De novo转录组分析 -- 分析笔记本

## 分析概览

本笔记本提供无参考基因组物种的转录组组装与分析完整流程，涵盖组装、质量评估、功能注释、编码区预测、表达定量和差异分析。

---

## 1. 转录组组装

### 1.1 Trinity组装 (推荐)

```bash
# 标准Trinity组装
Trinity --seqType fq \
  --left sample_R1.fq.gz \
  --right sample_R2.fq.gz \
  --max_memory 100G \
  --CPU 16 \
  --min_contig_length 200 \
  --output trinity_out_dir

# 多样本合并组装
Trinity --seqType fq \
  --samples_file samples.txt \
  --max_memory 200G \
  --CPU 32 \
  --output trinity_combined
```

**Trinity参数说明**：
- `--min_contig_length 200`: 最小contig长度，植物转录本建议200bp
- `--max_memory`: 最大内存，建议100-200G
- `--path_reinforcement_distance`: 路径强化距离，影响组装连续性

### 1.2 rnaSPAdes组装

```bash
# SPAdes转录组组装
rnaspades.py -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
  -t 16 -m 100 \
  -o rnaspades_out
```

### 1.3 混合组装策略

```bash
# 短读长 + 长读长混合组装
# Step 1: 短读长组装
Trinity --seqType fq --left R1.fq --right R2.fq --output short_asm

# Step 2: 长读长延长
LoRDEC -s short_asm/Trinity.fasta \
  -i long_reads.fq \
  -k 23 -s 3 \
  -o extended.fasta
```

---

## 2. 组装质量评估

### 2.1 BUSCO完整性评估

```bash
# 使用植物数据库评估
busco -i Trinity.fasta \
  -l embryophyta_odb10 \
  -m transcriptome \
  -o busco_output \
  -c 16

# 输出指标：
# Complete: 完整转录本比例
# Fragmented: 片段化转录本比例
# Missing: 缺失基因比例
```

**植物特异性考量**：
- 使用embryophyta_odb10数据库
- 完整度>80%为良好，>90%为优秀
- 多倍体物种可能有较高的重复基因BUSCO

### 2.2 TransRate质量评估

```bash
# TransRate综合评估
transrate --assembly Trinity.fasta \
  --left R1.fq --right R2.fq \
  --threads 16 \
  --output transrate_out

# 关键指标：
# n_seqs: 转录本数量
# gc_content: GC含量
# ortholog hit ratio: 直系同源命中率
```

### 2.3 N50和长度统计

```bash
# 使用TrinityStats.pl
TrinityStats.pl Trinity.fasta

# 辸出：
# Total trinity transcripts: 总转录本数
# Total assembled bases: 总碱基数
# Contig N50: N50值
# Median contig length: 中位数长度
```

---

## 3. 冗余去除

### 3.1 CD-HIT-EST去冗余

```bash
# 标准去冗余
cd-hit-est -i Trinity.fasta \
  -o Trinity_nr.fasta \
  -c 0.95 \
  -n 10 \
  -T 16 \
  -M 50000

# 参数说明：
# -c 0.95: 95%相似度聚类
# -n 10: 字符串长度（对应95%相似度）
```

### 3.2 Corset表达聚类

```bash
# 基于表达模式聚类
corset -i sample_counts.txt \
  -g sample_groups.txt \
  -n 0.5 \
  -o corset_out

# 输出cluster代表性转录本
```

### 3.3 EvidentialGene

```bash
# 综合去冗余流程
tr2aacds.pl -mrna Trinity.fasta

# 输出：
# okayalt.fa: 主要转录本
# doubtful.fa: 可疑转录本
```

---

## 4. 功能注释

### 4.1 Trinotate流程

```bash
# Step 1: ORF预测
TransDecoder.LongOrfs -t Trinity.fasta

# Step 2: BLAST搜索
blastp -query longest_orfs.pep \
  -db swissprot \
  -max_target_seqs 1 \
  -outfmt 6 \
  -evalue 1e-5 \
  -out blastp.out

# Step 3: Pfam搜索
hmmscan --cpu 8 \
  --domtblout pfam.domtblout \
  Pfam-A.hmm \
  longest_orfs.pep

# Step 4: 生成注释报告
Trinotate Trinity.fasta \
  --gene_trans_map Trinity.fasta.gene_trans_map \
  --transdecoder_pep longest_orfs.pep \
  --pfam_domain pfam.domtblout \
  --blastp blastp.out \
  > trinotate_annotation.tsv
```

### 4.2 植物特异性注释

```bash
# 使用植物数据库
# TAIR注释
blastp -query proteins.fa -db TAIR10_pep -out tair_blast.out

# PlantTFDB转录因子
hmmscan --tblout tf_domains.tbl PlantTFDB.hmm proteins.fa

# KEGG植物通路
ghostKOALA或KAAS在线注释
```

---

## 5. 编码区预测

### 5.1 TransDecoder

```bash
# 完整流程
TransDecoder.LongOrfs -t Trinity.fasta
TransDecoder.Predict -t Trinity.fasta

# 同源辅助预测
blastp -query longest_orfs.pep -db swissprot -out homology.out
TransDecoder.Predict -t Trinity.fasta --retain_blastp_hits homology.out
```

### 5.2 ORF质量过滤

```bash
# 过滤短ORF
awk 'length($2) >= 100' orfs.txt > filtered_orfs.txt

# 保留完整ORF
grep -E "complete|partial" orfs.txt
```

---

## 6. 表达定量

### 6.1 Salmon快速定量

```bash
# 构建索引
salmon index -t Trinity.fasta -i trinity_index --type quasi -k 31

# 定量
salmon quant -i trinity_index \
  -l A \
  -1 sample_R1.fq.gz \
  -2 sample_R2.fq.gz \
  -p 8 \
  -o sample_quant

# 合并多样本
salmon quantify -i trinity_index -l A -1 R1.fq -2 R2.fq -o quant1
# 使用tximport导入R
```

### 6.2 RSEM定量

```bash
# 构建参考
rsem-prepare-reference Trinity.fasta rsem_ref

# 定量
rsem-calculate-expression --paired-end \
  --num-threads 8 \
  sample_R1.fq sample_R2.fq \
  rsem_ref \
  sample_name
```

---

## 7. 差异表达分析

### 7.1 DESeq2分析

```r
library(tximport)
library(DESeq2)

# 导入Salmon定量结果
txi <- tximport(files, type="salmon", txOut=TRUE)

# 创建DESeq2对象
dds <- DESeqDataSetFromTximport(txi, colData, design=~condition)

# 差异分析
dds <- DESeq(dds)
res <- results(dds, lfcThreshold=1, alpha=0.05)

# 植物小样本策略：lfcThreshold=1
```

### 7.2 edgeR分析 (两重复)

```r
library(edgeR)

# 创建DGEList
y <- DGEList(counts=counts, group=group)

# 标准化
y <- calcNormFactors(y)

# 差异分析
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
```

---

## 8. 长非编码RNA分析

### 8.1 编码潜力评估

```bash
# CPC2评估
python CPC2.py -i transcripts.fa -o cpc2_output

# CNCI植物模型
CNCI.py -i transcripts.fa -m pl -o cnci_output

# FEELnc
FEELnc_filter.pl -i transcripts.fa --monoex=-1
```

### 8.2 lncRNA筛选标准

```bash
# 综合筛选
# 1. 长度 >= 200bp
# 2. 外显子数 >= 2 (可选)
# 3. CPC2 score < 0
# 4. CNCI non-coding
# 5. 无已知蛋白结构域
```

---

## 常见问题

### Q: 组装N50较低怎么办？
A: 增加测序深度，检查数据质量，考虑使用混合组装策略。

### Q: 转录本数量过多怎么办？
A: 使用更严格的去冗余参数，或使用Corset基于表达聚类。

### Q: 功能注释率低怎么办？
A: 使用近缘物种蛋白数据库，添加eggNOG-mapper注释。

### Q: 差异分析结果验证？
A: 使用RT-qPCR验证关键基因，检查生物学通路是否合理。

---

## 植物特异性考量

### 高杂合度物种
- 可能需要更高的测序深度
- 考虑单倍型分离组装

### 多倍体物种
- 同源基因区分困难
- 可能需要亚基因组特异性分析

### 非模式物种
- 功能注释依赖近缘物种
- 考虑使用通用数据库补充

---

## 参考

- Trinity: https://github.com/trinityrnaseq/trinityrnaseq
- TransDecoder: https://github.com/TransDecoder/TransDecoder
- Trinotate: https://github.com/Trinotate/Trinotate.github.io
- Salmon: https://combine-lab.github.io/salmon/
- BUSCO: https://busco.ezlab.org/
