# 植物 ChIP-seq 分析 — 特殊注意事项

## 植物 ChIP-seq 实验的特殊挑战

### 1. 植物材料的特殊性

#### 细胞壁与次生代谢物
- **细胞壁**: 植物细胞壁由纤维素、半纤维素、果胶和木质素组成，需要特殊裂解条件
- **次生代谢物**: 多酚、多糖等次生代谢物会干扰 ChIP 实验
  - 交联前使用真空浸润去除空气
  - 使用 PVPP (聚乙烯聚吡咯烷酮) 去除多酚
  - 在裂解液中添加 protease inhibitor cocktail
- **叶绿体**: 绿色组织中叶绿体 DNA 含量高，可能增加背景噪音

#### 组织类型
- **叶片**: 最常用的组织，叶绿体 DNA 污染需注意
- **根**: 获得足够生物量较困难
- **花/花序**: 发育阶段特异性，重复性可能受影响
- **愈伤组织**: 均一性较好，但可能与植物体内状态不同
- **幼苗**: 整株幼苗常用于组蛋白修饰研究

### 2. 基因组特性

#### 基因组大小与复杂度

| 物种 | 基因组大小 | 倍性 | ChIP-seq 测序深度建议 |
|------|-----------|------|---------------------|
| 拟南芥 | ~135 Mb | 二倍体 | 20-30M reads (TF), 40-50M (histone) |
| 水稻 | ~430 Mb | 二倍体 | 30-40M reads (TF), 50-70M (histone) |
| 玉米 | ~2.4 Gb | 古四倍体 | 50-80M reads (TF), 100M+ (histone) |
| 大豆 | ~1.1 Gb | 古四倍体 | 40-60M reads (TF), 80M+ (histone) |
| 小麦 | ~17 Gb | 六倍体 | >>100M reads |

#### 重复序列与多倍体
- 植物基因组中 **转座子 (TE)** 含量高 (玉米 > 85%)
- **多倍体**基因组中同源基因区分困难
- 使用 `-k 1` 限制唯一比对时可能丢失同源区域的信息
- 对多倍体物种，考虑亚基因组特异性分析

#### 基因组注释质量
- 模式物种注释较完善 (拟南芥、水稻)
- 非模式物种注释可能不完整
- GFF/GTF 文件质量影响 peak 注释准确性

### 3. 实验设计

#### 对照选择
- **Input DNA**: 去除超声打断后的总 DNA，标准对照
- **Mock IP** / **IgG**: 非特异性结合对照
- **植物特有**: 使用非转基因野生型作为转基因株系的对照

#### 生物学重复
- 推荐至少 **2 个生物学重复**
- 植物生长条件 (光照、温度、湿度) 需要严格控制
- 不同批次的植物材料可能导致批次效应

#### 抗体选择
- 植物组蛋白大部分保守，可使用商业抗体
  - H3K4me3, H3K27me3, H3K9ac, H3K27ac 有商业抗体可用
- 植物特异性组蛋白修饰 (如 H3K9me2) 需要验证抗体特异性
- 转录因子 ChIP 通常需要特异性抗体或转基因标签 (GFP/FLAG/HA)

### 4. 数据分析特殊考虑

#### Peak Calling
- **MACS2**: ChIP-seq peak calling 的金标准工具，支持 narrow peak (TF) 和 broad peak (组蛋白修饰) 模式
- **基因组大小参数**: 使用有效基因组大小 (mappable genome size)，不是总基因组大小
- **Input 归一化**: 对高重复序列基因组尤为重要
- **IR (Irreproducible Discovery Rate)**: IDR 框架用于评估重复一致性

#### 差异结合分析
- **归一化方法**: 考虑使用 spike-in (如 Drosophila chromatin) 进行归一化
- **多倍体 reads 分配**: 使用 `featureCounts` 的 multi-mapping 处理选项
- **Consensus peaks**: 使用多个重复的交集或至少 2/3 的 overlapping peaks

#### Motif 分析
- **植物 TF 数据库**:
  - PlantTFDB: http://planttfdb.gao-lab.org/
  - JASPAR Plants: https://jaspar.genereg.net/
  - CIS-BP: http://cisbp.ccbr.utoronto.ca/
  - PlantPAN: http://plantpan.itps.ncku.edu.tw/

### 5. 植物特有组蛋白修饰

| 修饰 | 功能 | Peak 类型 | 分布特征 |
|------|------|----------|---------|
| H3K4me3 | 活跃启动子 | Narrow | 基因 5‘ 端富集 |
| H3K27ac | 活跃增强子 | Broad | 基因间区和启动子 |
| H3K27me3 | 转录抑制 | Broad (宽域) | 广泛分布于基因体 |
| H3K9me2 | 转座子沉默 | Broad | TE 区域富集 (植物特有) |
| H3K36me3 | 转录延伸 | Broad | 基因 3’ 端富集 |
| H3K4me1 | 增强子标记 | Broad | 基因间区 |

### 6. 常用植物 ChIP-seq 资源

| 资源 | 描述 | 链接 |
|------|------|------|
| Plant Chromatin State Database | 植物染色质状态图谱 | http://systemsbiology.cau.edu.cn/chromstates/ |
| PlantDHS | 植物 DNase I 超敏位点 | http://plantdhs.org/ |
| ChIP-Hub | 植物 ChIP-seq 数据集成 | http://www.chiphub.org/ |
| PCSD | Plant Chromatin State Database | 多种植物表观基因组 |

### 7. 常见陷阱与解决方案

| 问题 | 原因 | 解决方案 |
|------|------|---------|
| 高背景信号 | 叶绿体 DNA 污染 | 比对后移除叶绿体 reads |
| Peak 数过少 | 抗体效率低 | 检查抗体特异性，优化 IP 条件 |
| 同源区域信号分散 | 多倍体/重复序列 | 允许 multi-mapping 或使用亚基因组特异性分析 |
| 组蛋白修饰边界不清 | Broad 域特征 | 使用 --broad 模式，调整 --broad-cutoff |
| 批次效应 | 不同批次的实验 | 使用 ComBat-seq 或 RUVseq 进行批次校正 |
