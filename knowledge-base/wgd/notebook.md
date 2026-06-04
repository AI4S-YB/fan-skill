# 全基因组复制分析 -- 分析笔记本

## 分析概览

本笔记本提供植物全基因组复制(WGD)分析的完整流程，涵盖KS分布分析、共线性分析、WGD检测、年龄估计和基因保留分析。

---

## 1. KS分布分析

### 1.1 WGDI KS分析

```bash
# WGDI配置文件
# wgdi.conf
[species]
species_name = target_species
genome = genome.fa
protein = protein.fa
gff = genes.gff
lens = chromosome.length

# 运行KS分析
wgdi -ks wgdi.conf

# 输出KS分布图
# ks_distribution.pdf
```

### 1.2 手动KS分析流程

```bash
# Step 1: 同源基因对鉴定
blastp -query proteins.fa -db proteins.fa \
  -evalue 1e-5 -outfmt 6 -out self_blast.out

# Step 2: 提取最佳匹配
awk '{if($1!=$2) print $1"\t"$2}' self_blast.out | \
  sort -k1,1 -k12,12g | \
  awk '!seen[$1]++' > best_hits.txt

# Step 3: 提取序列
python extract_pairs.py best_hits.txt proteins.fa pairs.fasta

# Step 4: 多序列比对
mafft --auto pairs.fasta > aligned.fasta

# Step 5: PAML yn00计算KS
yn00 < control_file

# Step 6: 绘制KS分布
Rscript plot_ks.R ks_values.txt
```

### 1.3 KS分布解读

| KS峰位置 | 大致年龄 | 解释 |
|---------|---------|------|
| 0-0.2 | <5 Mya | 近期WGD或种内复制 |
| 0.2-0.5 | 5-20 Mya | 较近期WGD |
| 0.5-1.0 | 20-50 Mya | 中等年龄WGD |
| 1.0-2.0 | 50-100 Mya | 较古老WGD |
| >2.0 | >100 Mya | 古老WGD |

**注意**: 精确年龄需要分子钟校准

---

## 2. 共线性分析

### 2.1 MCScanX分析

```bash
# 准备输入文件
# 基因位置文件格式: chr gene_id start end
awk '$3=="gene"{print $1"\t"$9"\t"$4"\t"$5}' genes.gff | \
  sed 's/ID=//;s/;.*//' > gene_location.txt

# 运行BLASTP
blastp -query proteins.fa -db proteins.fa \
  -evalue 1e-5 -outfmt 6 -out self.blast

# 运行MCScanX
MCScanx gene_location

# 查看结果
head gene_location.collinearity
```

**MCScanX输出解读**:
- `.collinearity`: 共线性区块列表
- `.tandem`: 串联重复基因
- `.wgd`: WGD来源的基因对
- `.segmental`: 片段复制基因对

### 2.2 WGDI共线性点图

```bash
# WGDI点图配置
wgdi -d wgdi.conf

# 输出点图直观展示WGD模式
# 点图中平行对角线表示WGD
```

### 2.3 共线性深度分析

```bash
# 统计共线性区块的拷贝数分布
python synteny_depth.py gene_location.collinearity

# 输出：
# 1:1 共线性区块数量
# 2:2 共线性区块数量 (暗示WGD)
# 4:4 共线性区块数量 (暗示两次WGD)
```

---

## 3. WGD检测与确认

### 3.1 多方法交叉验证

```bash
# 方法1: KS分布
# 检查KS分布是否有明显峰

# 方法2: 共线性分析
# 检查是否存在大规模共线性区块

# 方法3: 同源性深度
# 检查基因组内多拷贝同源基因比例

# 综合判断
python integrate_wgd_evidence.py ks_results.txt collinearity.txt depth_results.txt
```

### 3.2 区分共享WGD和独立WGD

```bash
# 比较多个物种的KS分布
# 如果两个物种的KS峰位置相近，可能是共享WGD

R
ks1 <- read.table("species1_ks.txt")
ks2 <- read.table("species2_ks.txt")
plot(density(ks1$V1), col="red")
lines(density(ks2$V1), col="blue")
```

---

## 4. WGD年龄估计

### 4.1 基于KS的年龄估计

```bash
# 使用分子钟校准
# 植物常用替代率: 1.5e-8 to 6e-9 substitutions/site/year

# 公式: Age = KS / (2 * substitution_rate)
python estimate_wgd_age.py --ks_peak 0.8 --rate 1.5e-8

# 输出估计年龄（百万年）
```

### 4.2 跨物种校准

```bash
# 结合分化时间校准
# 如果WGD峰在物种分化之前，则WGD是共享的

# 需要物种树和分化时间
python cross_species_wgd_dating.py \
  --species_tree tree.nwk \
  --ks_peaks ks_peaks.txt \
  --divergence_times divergence.txt
```

---

## 5. 祖先核型重建

### 5.1 多物种重建

```bash
# 使用ANCESTRALCHROM或类似工具
# 需要多个近缘物种的共线性数据

ancestralchrom \
  -i synteny_blocks.txt \
  -t species_tree.nwk \
  -o ancestral_output/

# 输出祖先染色体数和结构
```

### 5.2 单物种推断

```bash
# 基于基因组内共线性推断
# 识别嵌套共线性区块，推断祖先染色体

python infer_ancestral_karyotype.py \
  --collinearity gene_location.collinearity \
  --output ancestral_karyotype.txt
```

---

## 6. 基因保留分析

### 6.1 WGD基因对识别

```bash
# 从MCScanX结果提取WGD基因对
grep "WGD" gene_location.collinearity > wgd_pairs.txt

# 提取保留基因列表
awk '{print $2"\n"$3}' wgd_pairs.txt | sort -u > retained_genes.txt
```

### 6.2 保留基因功能分析

```bash
# GO富集分析
# 使用eggNOG-mapper注释
emapper.py -i retained_genes.fa -o annotation

# GO富集
python go_enrichment.py retained_genes.txt all_genes.txt go_annotation.txt
```

### 6.3 基因类别保留率

```python
# 计算不同功能类别的保留率
import pandas as pd

# 基因功能分类
categories = {
    'transcription_factor': [...],
    'kinase': [...],
    'photosynthesis': [...],
    # ...
}

# 计算每个类别的保留率
for cat, genes in categories.items():
    retained = len(set(genes) & set(retained_genes))
    total = len(genes)
    rate = retained / total
    print(f"{cat}: {rate:.2%}")
```

---

## 7. WGD后进化分析

### 7.1 基因丢失模式

```bash
# 分析WGD后的基因丢失
# 比较预期保留和实际保留

python gene_loss_analysis.py \
  --wgd_pairs wgd_pairs.txt \
  --current_genes current_genes.txt \
  --output loss_analysis.txt
```

### 7.2 表达分化分析

```bash
# 如果有RNA-seq数据
# 分析WGD基因对的表达差异

# 提取WGD基因对表达量
python extract_wgd_expression.py wgd_pairs.txt expression_matrix.txt

# 计算表达相关性
Rscript expression_correlation.R wgd_expression.txt
```

---

## 8. 植物特异性分析

### 8.1 古多倍体识别

```bash
# 一些"二倍体"实际是古多倍体
# 检查特征：
# 1. KS分布有古老峰
# 2. 大量基因有多拷贝
# 3. 复杂的共线性模式

python ancient_polyploid_detection.py ks.txt collinearity.txt
```

### 8.2 多次WGD历史

```bash
# 植物常有多次WGD
# 分析多次WGD的策略

# 1. 识别所有KS峰
Rscript identify_ks_peaks.R ks.txt

# 2. 区分不同时期的共线性区块
python stratify_synteny_blocks.py collinearity.txt ks_peaks.txt
```

---

## 常见问题

### Q: KS分布无明显峰？
A: 可能没有近期WGD，或组装/注释质量有问题。检查共线性分析结果。

### Q: KS峰太宽？
A: 可能代表多次WGD或连续的基因复制。尝试用混合模型分解。

### Q: 共线性区块太少？
A: 检查组装质量。高度碎片化的组装无法检测共线性。

### Q: 如何区分WGD和片段复制？
A: WGD产生全基因组范围的共线性，片段复制局限于局部区域。

---

## 参考

- WGDI: https://github.com/liu-shuijiu/WGDI
- MCScanX: https://github.com/wyp1125/MCScanX
- PAML: http://abacus.gene.ucl.ac.uk/software/paml.html
- ANCESTRALCHROM: https://github.com/JLSteenwyk/ANCESTRALCHROM
- CoGe SynFind: https://genomevolution.org/coge/SynFind.pl
