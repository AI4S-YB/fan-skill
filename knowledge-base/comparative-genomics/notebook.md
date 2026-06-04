# 比较基因组学分析 -- 分析笔记本

## 分析概览

本笔记本提供植物比较基因组学的完整分析流程，涵盖同源基因鉴定、基因组比对、共线性分析、WGD检测和进化推断。

---

## 1. 同源基因鉴定

### 1.1 OrthoFinder (推荐)

```bash
# 准备输入文件：每个物种的蛋白序列文件
# species1.fa, species2.fa, species3.fa ...

# 运行OrthoFinder
orthofinder -f protein_sequences/ \
  -t 16 \
  -a 8 \
  -S diamond

# 输出目录结构：
# Orthogroups/ - 同源基因组
# Orthogroup_Sequences/ - 每个组的序列
# Gene_Trees/ - 基因树
# Species_Tree/ - 物种树
```

**OrthoFinder参数说明**：
- `-f`: 输入目录
- `-t`: 线程数
- `-a`: 比对线程数
- `-S`: 搜索方法（diamond/blast）

### 1.2 结果解读

```bash
# 查看同源基因组统计
cat Orthogroups/Orthogroups.tsv

# 统计每个物种的直系同源基因数量
python ortholog_stats.py Orthogroups/Orthogroups.tsv
```

---

## 2. 基因组比对

### 2.1 minimap2

```bash
# 近缘物种比对
minimap2 -cx asm5 \
  reference.fa query.fa \
  > alignment.paf

# 更远缘物种
minimap2 -cx asm10 \
  reference.fa query.fa \
  > alignment.paf

# 转换为链格式
paftools.js call alignment.paf > alignment.chain
```

### 2.2 MUMmer

```bash
# NUCmer比对
nucmer --maxmatch -l 100 -c 500 \
  -p nucmer_output \
  reference.fa query.fa

# 生成delta文件
show-coords -rcl nucmer_output.delta > nucmer_output.coords

# 可视化
mummerplot --png nucmer_output.delta
```

### 2.3 LASTZ

```bash
# 高分歧物种比对
lastz reference.fa[multiple] query.fa \
  --output=alignment.maf \
  --format=maf
```

---

## 3. 共线性分析

### 3.1 MCScanX

```bash
# Step 1: 准备输入文件
# 基因位置文件：chr gene_id start end
awk '$3=="gene"{print $1"\t"$9"\t"$4"\t"$5}' annotation.gff3 | \
  sed 's/ID=//;s/;.*//' > gene_location.txt

# Step 2: 运行BLASTP
blastp -query proteins.fa -db proteins.fa \
  -evalue 1e-5 -outfmt 6 -out blast.out

# Step 3: 运行MCScanX
MCScanx gene_location

# 输出文件：
# gene_location.collinearity - 共线性区块
# gene_location.tandem - 串联重复
# gene_location.wgd - WGD基因对
```

### 3.2 JCVI

```python
# 使用JCVI进行共线性分析
from jcvi.compara.catalog import Ortholog
from jcvi.compara.synteny import AnchorFile

# 运行共线性分析
python -m jcvi.compara.catalog ortholog species1 species2

# 可视化
python -m jcvi.graphics.karyotype
```

### 3.3 Circos可视化

```bash
# 准备Circos输入文件
# 链接文件：chr1 start1 end1 chr2 start2 end2

# 生成Circos配置
circos -config circos.conf -outputdir circos_output/
```

---

## 4. WGD检测

### 4.1 WGDI

```bash
# WGDI完整流程
# Step 1: 准备BLASTP结果
blastp -query proteins.fa -db proteins.fa -out blast.out

# Step 2: 运行WGDI
wgdi -d conf_file  # 点图
wgdi -c conf_file  # 共线性
wgdi -k conf_file  # KS分布

# 参数配置文件示例
[species1]
species_name = species1
genome = species1.fa
protein = species1.pep
gff = species1.gff
```

### 4.2 KS分布分析

```bash
# 计算KS值
# 使用PAML yn00
yn00 < control_file

# 或使用KaKs_Calculator
KaKs_Calculator -i gene_pairs.axt -o kaks.txt -m MLWL

# 绘制KS分布
Rscript plot_ks_distribution.R kaks.txt
```

**KS分布解读**：
- 近期WGD：KS峰在0.1-0.3
- 中等年龄WGD：KS峰在0.5-1.0
- 古老WGD：KS峰在1.5-3.0

---

## 5. 染色体进化分析

### 5.1 祖先核型重建

```bash
# 使用ANCESTRALCHROM
ancestralchrom -i synteny_blocks.txt \
  -t species_tree.nwk \
  -o ancestral_output/

# 分析染色体重排
python analyze_rearrangements.py ancestral_output/
```

### 5.2 染色体数目进化

```bash
# 统计每个物种的染色体数目
for fa in *.fa; do
  echo -n "$fa: "
  grep "^>" $fa | wc -l
done

# 结合系统发育分析染色体数目变化
```

---

## 6. 基因家族进化

### 6.1 CAFE分析

```bash
# 准备输入文件
# 基因家族计数矩阵
# 物种树

# 运行CAFE
cafe5 -i gene_families.txt -t species_tree.nwk

# 输出：
# 显著扩张/收缩的基因家族
# 各分支的基因家族变化
```

---

## 7. 分歧时间估算

### 7.1 MCMCtree 贝叶斯分子钟分析

```bash
# MCMCtree是PAML包中的贝叶斯分子钟工具

# Step 1: 准备输入文件
# - 密码子/蛋白多序列比对 (phylip格式)
# - 带有标定点的物种树 (newick格式)

# 物种树示例（带化石标定点）
cat > species_tree_calibrated.nwk << EOF
(((species1:0.1,species2:0.1)'<0.5':0.2,species3:0.3):0.1,outgroup:0.5);
EOF

# Step 2: 配置mcmctree.ctl文件
cat > mcmctree.ctl << EOF
seed = -1
seqfile = alignment.phy
treefile = species_tree_calibrated.nwk
outfile = mcmctree_out
ndata = 1
seqtype = 2  # 0:codons, 1:amino acids, 2: nucleotides
usedata = 3  # 3: approximate likelihood calculation
clock = 2    # 2: independent rates (IR) model
RootAge = '<1.0'  # 根节点最大年龄约束
model = 0    # 0:JC69, 4:HKY85
alpha = 0.5  # Gamma distribution shape parameter
ncatG = 5    # Number of categories for Gamma
burnin = 100000
sampfreq = 10
nsample = 20000
EOF

# 运行MCMCtree
mcmctree mcmctree.ctl

# 输出文件：
# mcmctree_out.txt - 各节点分歧时间估计
# FigTree.tre - 可视化用树文件（含时间信息）
```

### 7.2 标定点选择指南

```bash
# 常用植物化石标定点来源：
# 1. 被子植物冠群年龄：约140-180 Mya
# 2. 单子叶植物冠群年龄：约110-130 Mya
# 3. 蔷薇类冠群年龄：约100-120 Mya
# 4. 物种特异性化石记录

# 标定点设置格式 (用于MCMCtree):
# '>0.15<0.25' 表示区间约束 0.15-0.25
# 'B(0.1,0.2)' 表示软边界 0.1-0.2
# 'L(0.05)' 表示最小年龄约束
```

**MCMCtree结果解读**：
- 各节点的后验平均年龄及95% HPD区间
- 检查ESS值确保MCMC收敛 (ESS > 200)
- 使用FigTree可视化带时间轴的树

---

## 8. 成对Ka/Ks选择压力分析

### 8.1 ParaAT + KaKs_Calculator 批量分析

```bash
# 准备输入文件
# 为每对物种准备：
# - 同源基因对列表（两列：物种A基因ID和物种B基因ID）
# - 蛋白序列文件
# - CDS序列文件

# Step 1: 提取直系同源基因对
cat orthogroups.tsv | while read line; do
  # 提取物种A和物种B的基因对
  echo "$line" | awk '{print $1"\t"$2}' >> gene_pairs.txt
done

# Step 2: 运行ParaAT
ParaAT.pl -h speciesA.pep -h speciesB.pep \
  -n speciesA.cds -n speciesB.cds \
  -p proc \
  -m mafft \
  -f axt \
  -k KaKs_Calculator \
  -o paraat_output/

# Step 3: 汇总结果
cat paraat_output/*.kaks > all_pairs.kaks

# Step 4: 检查结果
awk '{print $1"\t"$3/$4}' all_pairs.kaks | \
  sort -t$'\t' -k2 -nr | head -20
# 输出Ka/Ks最高的基因对（候选正选择基因）
```

### 8.2 手动成对Ka/Ks计算流程

```bash
# 单对基因的手动流程
# Step 1: 蛋白序列比对
mafft --auto gene_pair.pep > gene_pair_aligned.pep

# Step 2: 回译为密码子比对
pal2nal.pl gene_pair_aligned.pep gene_pair.cds \
  -output paml -nogap > codon_aligned.nuc

# Step 3: 计算Ka/Ks
KaKs_Calculator -i codon_aligned.nuc \
  -o kaks_output.txt \
  -m YN  # Yang-Nielsen方法

# 不同方法的选择：
# - YN: 适用于中等分歧物种
# - MA: 最大似然法，更精确但较慢
# - MLWL: M-L方法变体
# - GY: Goldman-Yang方法

# Step 4: 筛选显著正选择基因
awk '$5 > 1 && $6 < 0.05' kaks_output.txt > positively_selected.txt
# 条件：Ka/Ks > 1 且 P-value < 0.05
```

**Ka/Ks分析要点**：
- 优先使用高质量一对一同源基因
- 去除比对质量差的基因对
- Ka/Ks = 0.5-1.0 提示松弛选择
- 多个基因Ka/Ks > 1 可能指示适应性进化

---

## 9. LTR反转录转座子分析

### 9.1 LTRharvest de novo鉴定

```bash
# 安装GenomeTools
# conda install -c bioconda genomethools

# Step 1: 基因组索引
gt suffixerator -db genome.fasta \
  -indexname genome_index \
  -tis -suf -lcp -des -ssp -sds -dna

# Step 2: 运行LTRharvest
gt ltrharvest -index genome_index \
  -minlenltr 100 \
  -maxlenltr 7000 \
  -mindistltr 1000 \
  -maxdistltr 25000 \
  -similar 80 \
  -mintsd 4 -maxtsd 20 \
  -vic 10 \
  -seed 20 \
  -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 \
  -out ltrharvest_output.fasta \
  -outinner ltrharvest_inner.fasta \
  > ltrharvest_results.gff3

# 参数说明：
# -minlenltr/maxlenltr: LTR最小/最大长度
# -mindistltr/maxdistltr: LTR对最小/最大距离
# -similar: LTR对相似度阈值 (%)
# -mintsd/maxtsd: TSD最小/最大长度
# -vic: LTR识别所需邻近匹配数
```

### 9.2 LTR_retriever 过滤与插入时间估计

```bash
# Step 1: 准备LTR_retriever输入
# 需要：基因组文件 + LTRharvest结果 + TEsorter输出

# Step 2: 运行TEsorter鉴定LTR类型
TEsorter -db rexdb -p 16 \
  ltrharvest_output.fasta \
  > ltrharvest_tesorter.tsv

# Step 3: 运行LTR_retriever
LTR_retriever -genome genome.fasta \
  -inharvest ltrharvest_results.gff3 \
  -infinder ltrfinder_results.txt \
  -threads 16 \
  -u ${species_substitution_rate}

# 不同物种替代率参考：
# - 拟南芥 (Arabidopsis): 7.0e-9
# - 水稻 (Oryza sativa): 1.3e-8
# - 玉米 (Zea mays): 3.0e-8
# - 大豆 (Glycine max): 6.0e-9

# Step 4: 提取完整LTR元件
cat genome.fasta.mod.pass.list | while read id; do
  samtools faidx genome.fasta $id >> intact_ltr.fasta
done

# Step 5: LTR插入时间分布分析
Rscript -e '
library(ggplot2)
ltr_age = read.table("genome.fasta.mod.pass.list.age", header=TRUE)
ggplot(ltr_age, aes(x=Insertion_Time_MYA)) +
  geom_histogram(bins=50, fill="steelblue") +
  xlab("LTR Insertion Time (MYA)") +
  ylab("Count") +
  ggtitle("LTR Insertion Time Distribution")
'
```

**LTR分析解读**：
- 插入时间峰反映LTR爆发事件
- 近期峰可能与物种形成或驯化相关
- 古老LTR通常在基因组中累积突变较多
- LTR占植物基因组比例可反映基因组大小演化

---

## 10. 正选择分析

### 10.1 PAML codeml

```bash
# 准备密码子比对和树文件

# 分支模型检测特定谱系正选择
# codeml.ctl配置文件
seqfile = codon_alignment.phy
treefile = species_tree.nwk
outfile = mlc
runmode = 0
model = 2  # branch model
NSsites = 0

codeml codeml.ctl

# 位点模型检测位点正选择
# 修改配置
model = 0
NSsites = 2  # sites model (M2 vs M1)
```

### 10.2 HyPhy

```bash
# BUSTED方法
hyphy busted --alignment codon_alignment.fa --tree species_tree.nwk

# aBSREL方法（分支特异性）
hyphy absrel --alignment codon_alignment.fa --tree species_tree.nwk

# FEL方法（位点特异性）
hyphy fel --alignment codon_alignment.fa --tree species_tree.nwk
```

---

## 11. 植物特异性分析

### 11.1 多倍化历史分析

```bash
# 植物基因组常有多次WGD
# 使用KS分布和共线性联合分析

# 区分共享WGD和独立WGD
python compare_wgd.py species1_ks.txt species2_ks.txt
```

### 11.2 核型进化

```bash
# 分析染色体融合/分裂
# 基于共线性区块的方向和位置

python analyze_fission_fusion.py synteny_blocks/
```

---

## 12. 结果整合

### 12.1 多分析整合

```bash
# 整合同源基因、共线性、WGD结果
# 生成进化历史报告

Rscript integrate_results.R \
  orthogroups.tsv \
  collinearity.txt \
  ks_distribution.txt \
  species_tree.nwk
```

---

## 常见问题

### Q: 同源基因鉴定结果太多怎么办？
A: 调整e-value阈值或使用更严格的一致性阈值。

### Q: 共线性区块太多？
A: 增加最小基因数阈值，减少噪音。

### Q: KS分布无明显峰？
A: 可能没有近期WGD，或组装/注释质量有问题。

### Q: 不同方法WGD推断不一致？
A: 综合多种证据，包括共线性、KS分布、基因树一致性。

---

## 参考

- OrthoFinder: https://github.com/davidemms/OrthoFinder
- MUMmer: https://mummer4.github.io/
- MCScanX: https://github.com/wyp1125/MCScanX
- WGDI: https://github.com/liu-shuijiu/WGDI
- CAFE: https://github.com/hahnlab/CAFE5
- PAML: http://abacus.gene.ucl.ac.uk/software/paml.html
- HyPhy: https://hyphy.org/
- Circos: http://circos.ca/
