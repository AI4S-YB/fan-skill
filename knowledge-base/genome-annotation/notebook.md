# 植物基因组注释：分析笔记

这份笔记帮助你像一个有经验的植物基因组学家一样思考和做决策。

## 在开始之前：理解你的基因组组装

`inspect_data.sh` 产出了 data_profile.json。在注释之前，先评估组装质量。

### 组装质量评估

查看 `assembly_n50` 和 `contig_count` 字段：

- **N50 > 1Mb**：染色体水平组装，适合完整注释。可以获得接近参考基因组的基因集合。
- **N50 100Kb–1Mb**：scaffold 水平组装，可以做基因预测，但某些基因可能被 scaffold 边界打断。
- **N50 < 100Kb**：片段化组装。基因预测的质量会受影响——预期更多不完整基因模型和融合基因。用 BUSCO 评估后谨慎解读。

### 这个基因组多大？

植物基因组大小变化巨大：

- **< 200Mb**：小基因组（如拟南芥 135Mb）。注释速度快，计算资源需求低。
- **200Mb–1Gb**：中等基因组（如水稻 430Mb，大豆 1.1Gb）。标准注释流程。
- **1Gb–5Gb**：大基因组（如玉米 2.4Gb，大麦 5Gb）。重复序列占比高（>70%），考虑分区注释或使用集群。
- **> 5Gb**：超大基因组（如小麦 17Gb，松树 22Gb）。必须分区注释。优先注释低拷贝区域。

### 重复序列含量

植物基因组以重复序列为主。LTR 逆转录转座子是最大组分。

- 重复序列 < 30% → 基因组较小或经历了重复序列清除
- 重复序列 30%–60% → 典型植物基因组水平
- 重复序列 > 60% → 常见于大基因组。RepeatMasker 是最耗时的步骤。

## 选择你的基因预测策略

### 有 RNA-seq 数据时的最佳方案

BRAKER3 是目前植物基因组注释的最优选择：

1. 用 HISAT2 将 RNA-seq reads 比对到屏蔽重复的基因组
2. 用 StringTie2 组装转录本
3. BRAKER3 利用 RNA-seq 比对和同源蛋白 hints 预测基因

RNA-seq 证据可以：
- 精确定位外显子-内含子边界
- 区分可变剪切体
- 验证基因表达的时空特异性

**覆盖度建议**：至少 50M paired-end reads，最好来自多种组织（叶、根、花、种子、胁迫处理）。

### 没有 RNA-seq 数据时

同源蛋白 hints 是退而求其次的选择：

1. 从 OrthoDB 下载植物蛋白集
2. 用 ProtHint 生成蛋白 hints
3. BRAKER3 利用蛋白 hints 预测基因

这种策略的局限：
- 只能找到保守基因
- 物种特有基因或快速进化基因可能漏掉
- 基因边界可能不够精确

### 多倍体基因组的注意事项

多倍体植物（小麦、棉花、油菜）面临特殊挑战：
- 亚基因组间的同源基因难以区分
- 对每个亚基因组单独注释，或使用亚基因组特异的 RNA-seq
- 注释完成后检查亚基因组的基因数是否大致平衡
- 高度相似的旁系同源基因可能导致基因融合

## 功能注释：不止是"跑个流程"

### EggNOG-mapper 能告诉你什么

- GO 术语：分子功能、生物学过程、细胞组分
- KEGG 通路：代谢通路映射
- COG 分类：功能类别概览
- 直系同源群归属

### InterProScan 的价值

InterProScan 提供 EggNOG 没有的信息：
- Pfam 结构域注释
- 跨膜区和信号肽预测
- 亚细胞定位线索
- 酶学委员会（EC）编号

**两套工具互补使用**是最佳实践。EggNOG 快但粒度粗，InterProScan 慢但信息丰富。

## 注释完整性评估

### BUSCO 的使用

BUSCO 蛋白质模式（`--mode proteins`）评估你预测的蛋白集合包含多少保守单拷贝直系同源基因。

解读 BUSCO 结果：
- **Complete > 90%**：优秀。注释质量接近参考基因组水平。
- **Complete 80%–90%**：良好。大部分保守基因已捕获。
- **Complete 70%–80%**：可接受。可能漏掉了一些基因。
- **Complete < 70%**：需要改进。考虑增加 RNA-seq 证据或调整 BRAKER3 参数。

### 与其他物种比较

不要孤立地看 BUSCO 分数。与近缘物种的已发表注释比较：
- 基因数量是否在预期范围内？
- 平均基因长度和外显子数量是否合理？
- 如果基因数量异常多（如 2 倍于近缘种），可能存在假阳性预测或基因组重复未被屏蔽

---

## 常见陷阱

### 基因数量异常多

可能原因：
1. 重复序列未充分屏蔽 → 转座子被误标为基因
2. 单倍型组装导致的冗余 → 同一位点在两个单倍型上各预测了一次
3. 可变剪切被计为独立基因 → 检查基因座（locus）数量而非仅仅转录本数量

### 基因数量异常少

可能原因：
1. BRAKER3 参数过于严格（如 evidence 数量阈值）
2. RNA-seq 覆盖度不足或组织来源单一
3. 基因组组装片段化导致基因被打断

### 基因间区太短

植物基因间区通常较短（平均 2-5Kb）。如果出现大量重叠基因：
- 检查是否是转座子被错误注释为基因
- 检查 UTR 区域是否被正确预测（BRAKER3 默认不预测 UTR）

### 非植物物种误用植物数据库

确保所有使用的数据库是植物特异版本：
- OrthoDB：使用 `Viridiplantae` 层级
- BUSCO：使用 `embryophyta_odb10` 或更特异的数据集
- InterProScan：不需要特殊配置，但结果解读应参考植物文献

---

## 基因预测工具选择：BRAKER3 vs GeMoMa vs MAKER

三个工具各有适用场景，选择取决于你拥有的证据类型。

### BRAKER3：有 RNA-seq 时的最优选择

**适用场景**：
- 拥有高质量 RNA-seq 数据（多组织、有重复）
- 需要自动化、可重复的流程
- 目标是获得包含可变剪切信息的完整基因模型

**优势**：自动训练 AUGUSTUS 和 GeneMark，无需人工干预；整合 RNA-seq 和蛋白 hints
**局限**：对 RNA-seq 质量敏感；训练过程不透明，调参空间有限

### GeMoMa：有近缘物种参考时的最佳选择

**适用场景**：
- 没有 RNA-seq 数据但近缘物种（同属或同科）有高质量注释
- 需要跨物种同源基因的保守结构
- 多倍体基因组（支持亚基因组特异的参考）

**优势**：基于保守基因结构，对外显子-内含子边界预测精确；可整合多个参考物种
**局限**：只能找到保守基因；物种特异性基因和快速进化基因会漏掉；需要近缘参考

### MAKER：需要最大灵活性和人工审编时

**适用场景**：
- 有异质证据来源（RNA-seq + 蛋白 + EST + 从头预测），需要手动整合
- 社区注释项目，需要人工审编
- 对注释质量要求极高（如作为参考基因组发表）

**优势**：整合所有证据类型；可人工审编和调整；注释结果可溯源
**局限**：配置复杂；需要大量人工干预；计算时间长；不适合高通量批处理

### 决策总结

| 证据 | 推荐工具 | 替代方案 |
|------|---------|---------|
| 高质量 RNA-seq | BRAKER3 | MAKER |
| 近缘物种参考蛋白 | GeMoMa | BRAKER3 + ProtHint |
| RNA-seq + 蛋白 | BRAKER3 + ProtHint | MAKER |
| 仅有从头预测 | GeMoMa + 远缘参考 | BRAKER3 + 远缘蛋白 |
| 多倍体基因组 | GeMoMa（分亚基因组） | BRAKER3 + 亚基因组 RNA-seq |

---

## RNA-seq 证据质量评估

不要盲目使用任何 RNA-seq 数据。在送入 BRAKER3 之前，评估以下指标。

### 链特异性

BRAKER3 默认假设 RNA-seq 数据是**非链特异性**的。如果你的数据是链特异性的（如 dUTP 方法），需要在比对时使用 `--rna-strandness` 参数，否则外显子-内含子边界的预测会出错。检查方法：使用 `infer_experiment.py`（RSeQC）推断链方向。

### 测序深度

- **最低要求**：每个组织 30M paired-end reads
- **推荐深度**：每个组织 50-80M paired-end reads
- **多组织合并**：至少涵盖叶、根、花/花序。胁迫处理（干旱、高温、病害）可额外增加稀有转录本的覆盖
- **警告**：深度不足（< 20M）会导致基因预测不完整，特别是低表达基因和长基因

### 组织覆盖度

植物基因表达具高度组织特异性。如果只有单一组织的 RNA-seq（如仅叶片），预期会漏掉：
- 根特异性次生代谢基因
- 花特异性转录因子
- 种子储存蛋白基因
- 胁迫响应基因

**最少推荐组织组合**：叶片 + 根 + 花/花序 + 一个胁迫处理（如有）。

### RNA 质量（RIN 值）

降解的 RNA 会导致 3' 偏向，影响基因模型完整性。RIN > 7 是可接受的阈值。对于植物组织（尤其是富含多糖多酚的根和果实），RIN 可能偏低（5-7），但即使如此，如果建库方案中有 polyA 捕获，仍可尝试使用。

---

## 重复序列屏蔽策略：植物基因组的关键

植物基因组中 30%-90% 是重复序列。如果不屏蔽，转座子会被误标为基因，导致基因数量大幅膨胀。

### 同源屏蔽 vs 从头屏蔽

| 策略 | 工具 | 适用场景 |
|------|------|---------|
| 同源屏蔽 | RepeatMasker + Repbase | 模式物种，重复库含该物种或近缘物种 |
| 从头屏蔽 | RepeatModeler2 + RepeatMasker | 非模式物种，无现成重复库 |
| 混合屏蔽 | RepeatModeler2 建库 + RepeatMasker + Repbase | 推荐方案：结合两者优势 |

**推荐流程**：
1. 使用 RepeatModeler2 对基因组进行从头重复序列预测，构建物种特异性重复库
2. 将物种特异库与 Repbase 植物部分合并
3. 使用 RepeatMasker 对基因组进行屏蔽
4. 为屏蔽后的基因组创建 BRAKER3 的软屏蔽索引（重复区域用小写字母标记）

### 物种特异性重复库的价值

- 非模式植物通常不在 Repbase 中有专门条目
- RepeatModeler2 从基因组本身学习重复序列，能发现物种特异的 TEs
- 对 LTR 逆转录转座子（植物最大的重复组分），使用 `LTRharvest` + `LTRdigest` 构建高质量 LTR 库
- 注意：RepeatModeler2 需要足够的计算资源（大基因组可能需要 > 100GB RAM）

### 软屏蔽 vs 硬屏蔽

- **软屏蔽**（推荐）：重复区域用小写字母标记，比对工具（HISAT2, minimap2）可以考虑这些区域，但不优先
- **硬屏蔽**：重复区域替换为 N，彻底排除
- BRAKER3 通常使用**软屏蔽**基因组，保留重复边缘区域可能含有的真实基因

---

## BUSCO 深入解读：分数如何影响下游使用

### 不同完整性分数的下游影响

| BUSCO Complete | 解读 | 下游使用建议 |
|---------------|------|------------|
| > 95% | 接近参考质量 | 适合发表、比较基因组学、基因家族分析 |
| 90-95% | 高质量 | 适合大多数分析，包括 GWAS 和 QTL 定位的候选基因注释 |
| 80-90% | 良好 | 可用于功能注释和通路分析，但进化分析需谨慎 |
| 70-80% | 可接受 | 通路分析基本可用，但"缺失基因"可能影响结论 |
| < 70% | 需改进 | 不要用于任何涉及基因集完整性的分析 |

### 不要只看 Complete 分数

- **Duplicated 比例**：二倍体物种 D < 10% 正常，多倍体中 D 升高是预期现象。但如果 D > 30%，可能不是多倍体而是组装冗余（两个 haplotypes 各自组装了）
- **Fragmented 比例**：F > 15% 提示基因模型不完整，可能是组装碎片化或 RNA-seq 覆盖不足导致的
- **Missing 比例**：M > 10% 提示有系统性缺失，检查 BUSCO lineage 选择是否合适（如非十字花科却用了 brassicales_odb10）

### 跨 lineage 比较

同时运行多个 BUSCO lineage（如同时跑 embryophyta_odb10 和物种所属目的 lineage）。通用 lineage 告诉你"缺了什么"，特化 lineage 告诉你"在近缘物种中缺了什么"。

---

## 非模式植物注释策略

大多数植物研究者不用拟南芥。当你的物种与任何已注释物种的亲缘关系 > 50Mya 时，需要调整策略。

### 跨物种蛋白证据的选择

使用 ProtHint 时，蛋白证据的进化距离影响预测质量：

- **同科**：保守度高，蛋白 hints 准确。优先使用。
- **同目**：可使用，但准确性下降。建议结合 RNA-seq（如有少量）。
- **同纲**：仅能捕获保守的看家基因。预期会漏掉 30-50% 的基因。
- **跨门**：不推荐。预测结果不可靠。

### 非模式植物 RNA-seq 策略

如果没有足够经费做多组织深度 RNA-seq：
1. **混合样本策略**：将多个组织的 RNA 等量混合后测序（降低成本但丢失组织特异性表达信息），加上一个深度测序的参考组织
2. **SRA 挖掘**：检查 NCBI SRA 中是否有同一属或科的 RNA-seq 数据可借用以提供额外 hints
3. **Iso-Seq 补充**：对非模式物种，少量 PacBio Iso-Seq（1-2 SMRT cells）可以提供全长转录本，大幅提高基因模型精度

### 从头训练基因预测器

如果没有任何近缘参考和 RNA-seq：
1. 使用 GeneMark-ES 进行从头基因预测（无需训练集）
2. 使用蛋白同源 hints（即使是远缘的）辅助
3. 结果应标记为 "预测基因" 而非 "注释基因"，报告应诚实说明局限性

---

## 非编码RNA (ncRNA) 注释

### ncRNA 注释的必要性

基因组注释不仅包括蛋白质编码基因，还包括各类非编码RNA基因。植物基因组中含有丰富的ncRNA：

- **rRNA**: 核糖体RNA (18S, 5.8S, 28S, 5S) — 在所有真核生物中高度保守
- **tRNA**: 转运RNA — 每个基因组通常有数百至数千个tRNA基因
- **snRNA/snoRNA**: 剪切体和核仁小RNA — 参与mRNA剪接和rRNA修饰
- **miRNA**: microRNA — 约21-24nt的调控RNA，在植物发育和胁迫响应中关键

### rRNA预测 (Barrnap)

```bash
# Barrnap 用于rRNA基因预测，支持真核生物rRNA鉴定
# 输出为GFF3格式
barrnap --kingdom euk --threads 8 genome.fasta > rrna.gff3

# 结果解读:
#   - 18S rRNA: 细胞质核糖体小亚基
#   - 5.8S rRNA: 细胞质核糖体大亚基 (与28S形成复合体)
#   - 28S rRNA: 细胞质核糖体大亚基
#   - 5S rRNA: 独立转录的核糖体大亚基组分
# 预期: 植物基因组中rRNA基因通常以串联重复阵列存在 (数百至数千拷贝)
```

### tRNA预测 (tRNAscan-SE)

```bash
# tRNAscan-SE 为tRNA基因预测的金标准工具
# 支持真核生物模式，输出GFF3格式
tRNAscan-SE -o trna.txt --gff3 trna.gff3 -E genome.fasta

# 参数说明:
#   -o: 文本格式的详细输出 (含二级结构预测)
#   --gff3: GFF3格式的基因组注释输出
#   -E: 真核生物模式 (搜索真核生物特异性tRNA特征)
# 
# 结果解读:
#   - 假基因: tRNAscan-SE 会标注可能的tRNA假基因
#   - 内含子: 植物tRNA基因可能含内含子 (tRNAscan-SE 会标记)
#   - 抑制tRNA: 标注可能的抑制tRNA (anticodon突变)
```

### 其他ncRNA鉴定 (Infernal cmscan + Rfam)

```bash
# Infernal cmscan 基于协方差模型搜索Rfam数据库
# 用于鉴定 snRNA, snoRNA, miRNA 等非编码RNA
cmscan --cut_ga --rfam --nohmmonly --tblout rfam.tblout \
  --fmt 2 --clanin Rfam.clanin Rfam.cm genome.fasta

# 参数说明:
#   --cut_ga: 使用Rfam curated的gathering threshold (高置信度)
#   --rfam: 使用Rfam加速选项 (fast mode)
#   --nohmmonly: 不运行HMM-only搜索 (仅使用CM模型)
#   --tblout: 表格格式输出文件
#   --fmt 2: 输出格式 (2 = 表格格式)
#   --clanin: Rfam clan信息文件 (用于去冗余)
#
# 注意: 此步骤需要Rfam数据库 (https://rfam.org/)
#   cmscan运行速度较慢 (大基因组可能需要数天)，建议在集群上运行
#   如不需要全面ncRNA注释，可仅运行Barrnap + tRNAscan-SE
```

### ncRNA注释的优先级

| ncRNA类型 | 工具 | 必需程度 | 计算时间 | 说明 |
|-----------|------|---------|---------|------|
| rRNA | Barrnap | 必需 | 分钟级 | 所有基因组注释都应包含rRNA |
| tRNA | tRNAscan-SE | 必需 | 分钟级 | tRNA是翻译机器的核心组分 |
| snRNA/snoRNA | Infernal cmscan | 推荐 | 小时-天 | 对理解剪接机制重要 |
| miRNA | Infernal cmscan 或专门工具 | 可选 | 天级 | 需要专门的miRNA预测流程更准确 |
| lncRNA | 专门流程 (CNCI+CPC2+FEELnc) | 可选 | 小时级 | 需要RNA-seq数据支持 |

---

## PASA基因结构更新与UTR鉴定

### PASA的价值

PASA (Program to Assemble Spliced Alignments) 是基因组注释精修的核心工具，主要用于：

1. **基因结构更新**: 基于RNA-seq转录本更新预测基因的外显子-内含子边界
2. **UTR鉴定**: 提供5'UTR和3'UTR的注释（AUGUSTUS/BRAKER3默认不预测UTR）
3. **可变剪接验证**: 通过与转录本比对验证和发现可变剪接事件
4. **注释比较**: 比较新旧注释版本的差异，追踪基因模型的改进

### PASA alignAssembly (转录本比对与基因结构更新)

```bash
# 将Trinity/StringTie转录本比对到基因组
# 用于更新基因结构注释、鉴定UTR和验证可变剪接
PASA/Launch_PASA_pipeline.pl \
  -c alignAssembly.config \
  -C \              # 创建MySQL数据库
  -R \              # 运行比对和组装
  -g genome.fasta \ # 基因组序列
  -t Trinity.fasta \ # 转录本FASTA (Trinity组装或StringTie输出)
  --ALIGNERS blat,gmap \  # 使用BLAT和GMAP两种比对工具
  --CPU 16

# alignAssembly.config 配置文件示例:
#   DATABASE=<mysql_database_name>
#   MYSQL_RW_USER=access
#   MYSQL_RW_PASSWORD=access
#   MYSQL_RW_HOST=localhost
#   validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
#   validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95
```

### PASA annotCompare (注释比较与UTR更新)

```bash
# 比较和更新已有注释的UTR边界
# 此步骤需要先完成alignAssembly
PASA/scripts/Load_Current_Gene_Annotations.dbi \
  -c annotCompare.config \  # 配置文件 (MySQL连接信息)
  -g genome.fasta \         # 基因组序列
  -P existing.gff3          # 需要更新的现有注释

# 运行比较
PASA/scripts/Launch_PASA_pipeline.pl \
  -c annotCompare.config \
  -A \  # 运行annotCompare模式
  -g genome.fasta \
  -t Trinity.fasta

# 输出:
#   - pasa_assembly.gff3: 更新后的基因注释 (含UTR)
#   - pasa_assembly.fasta: 更新后的转录本序列
#   - annotationCompare.report: 比较报告 (新旧注释差异)
```

### PASA使用注意事项

- **MySQL依赖**: PASA需要MySQL数据库存储中间结果。确保MySQL服务可用，用户有建库权限
- **转录本质量**: 输入的转录本应该是高质量的 (建议使用Trinity组装 + StringTie合并)
- **计算时间**: 大型转录本集 (> 100K 转录本) 的比对可能需要1-2天
- **内存需求**: PASA的MySQL实例在大型转录本集时可能需要 > 16GB RAM
- **与BRAKER3的配合**: BRAKER3输出不含UTR → PASA补充 → MAKER2整合含UTR的最终注释

---

## Diamond蛋白同源搜索补充

### 为什么需要独立的Diamond搜索

EggNOG-mapper虽然功能全面，但其内部使用MMseqs2进行同源搜索，数据库覆盖度有限。独立的Diamond搜索提供互补信息：

- **NR数据库**: NCBI非冗余蛋白数据库覆盖面最广 (含GenBank CDS翻译、RefSeq、Swiss-Prot、PRF、PDB)
- **Swiss-Prot数据库**: 高质量手工注释的蛋白数据库，注释信息最可靠
- **互检验证**: EggNOG和Diamond的top hit是否一致？不一致时可能存在功能分歧，需要深入分析

### NR数据库搜索

```bash
# Diamond blastp 搜索NR数据库
# 使用更敏感的blastp模式 (--sensitive)
diamond blastp \
  --db nr \
  -q final_proteins.faa \
  --evalue 1e-5 \
  --max-target-seqs 1 \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  --threads 16 \
  --sensitive \
  -o diamond_nr.out

# 输出格式说明 (outfmt 6 自定义列):
#   qseqid: 查询序列ID (你的蛋白)
#   sseqid: 目标序列ID (NR数据库条目)
#   pident: 相似度百分比
#   length: 比对长度
#   stitle: 目标序列的全称 (含物种和功能描述)
```

### Swiss-Prot数据库搜索

```bash
# Diamond blastp 搜索Swiss-Prot数据库
diamond blastp \
  --db swissprot \
  -q final_proteins.faa \
  --evalue 1e-5 \
  --max-target-seqs 1 \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  --threads 16 \
  --sensitive \
  -o diamond_swissprot.out

# Swiss-Prot注释特点:
#   - 所有条目经过人工审编 (curated)
#   - 功能注释经过实验验证
#   - 包含蛋白功能、亚细胞定位、翻译后修饰等信息
#   - 条目数量远小于NR (Swiss-Prot ~570K, NR > 500M)
```

### Diamond vs EggNOG-mapper 互补策略

```bash
# 建议的整合策略:
# 1. EggNOG-mapper → 获得 GO/KEGG/COG 系统分类
# 2. Diamond NR → 获得最广泛的同源匹配，用于物种分类和功能描述
# 3. Diamond Swiss-Prot → 获得最可靠的功能注释 (curated)
# 4. 三者结果整合: Swiss-Prot > NR > EggNOG 的优先级进行功能分配
```

### Diamond结果解读

- **相似度 (pident)**: > 70% 通常表示可靠的同源关系; 40-70% 需谨慎解读; < 40% 可能是不同蛋白家族
- **比对覆盖率**: 如果query coverage < 50%，可能是部分结构域匹配或多结构域蛋白的部分匹配
- **top hit物种**: 检查top hit的物种是否与你的物种亲缘关系合理
- **与EggNOG一致性**: 标记两种工具结果不一致的基因，它们可能是快速进化的物种特异性基因

---

## 额外常见陷阱

### 大基因组基因预测膨胀

玉米、大麦、小麦等大基因组（> 2Gb）中 > 80% 是重复序列。即使使用 RepeatMasker，仍可能：
- LTR 逆转录转座子的内部序列被误标为"逆转录酶基因"
- helitron 捕获的基因片段被预测为完整基因
- **对策**：对基因预测结果通过 Pfam 结构域过滤 — 如果预测的"基因"只含转座子相关结构域（如 integrase, reverse transcriptase, gag），丢弃

### 假基因过滤

植物基因组中假基因的处理：
- 与功能基因相比，假基因可能有：提前终止密码子、移码突变、不完整的结构域
- BRAKER3 可能将假基因注释为多个"碎裂基因"（fragmented genes）
- **过滤步骤**：注释完成后，用 `linclust` 或 `cd-hit` 聚类，标注与已知 TE 蛋白相似的序列，标记截断蛋白（< 50% 同源蛋白长度）

### 多倍体同源基因注释

多倍体植物中 homeolog 的注释是核心挑战：
- **亚基因组单独注释**：如果已区分亚基因组（如小麦 A/B/D），对每个亚基因组分别运行 BRAKER3，然后用同源关系链接 homeolog
- **亚基因组不平衡**：一个亚基因组注释到 40,000 个基因而另一个只有 30,000 个 → 检查是否有亚基因组偏倚的 RNA-seq 覆盖或蛋白证据
- **Homeolog 融合**：如果两个 homeolog 位于高相似性区域，BRAKER3 可能将它们合并为一个"超级基因"。检查异常长的基因（> 20Kb CDS），它们可能是融合产物
