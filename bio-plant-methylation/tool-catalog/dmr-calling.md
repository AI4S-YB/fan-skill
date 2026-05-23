# DMR 检测 -- 工具目录

## 概述

差异甲基化区域（DMR）检测用于鉴定不同条件之间（如处理 vs 对照、不同组织、发育阶段）甲基化水平显著变化的基因组区域。植物中需要在 CG、CHG 和 CHH 三种上下文分别进行 DMR 分析。

## 推荐工具

### 1. DSS (Dispersion Shrinkage for Sequencing) -- 金标准

**描述**: DSS 是专门为亚硫酸盐测序数据设计的 DMR 检测 R 包，基于贝叶斯分层模型和平滑化方法。不需要生物学重复也可以运行，但对有重复的实验统计效能更高。

**安装**:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DSS")
```

**完整分析流程**:

```r
library(DSS)

# 步骤 1: 读取甲基化数据
# 需要各样本的 cytosine report 格式文件
sample1 <- read.table("sample1_CG_report.txt",
  header = TRUE, sep = "\t")
sample2 <- read.table("sample2_CG_report.txt",
  header = TRUE, sep = "\t")
sample3 <- read.table("sample3_CG_report.txt",
  header = TRUE, sep = "\t")
sample4 <- read.table("sample4_CG_report.txt",
  header = TRUE, sep = "\t")

# 步骤 2: 构建 BSseq 对象
BSobj <- makeBSseqData(
  list(sample1, sample2, sample3, sample4),
  c("ctrl_rep1", "ctrl_rep2", "trt_rep1", "trt_rep2")
)

# 查看数据
BSobj
```

**参数说明**:
- `makeBSseqData()`: 将多个样本的甲基化数据合并为 BSseq 对象
- 样本名称列表与数据文件一一对应

```r
# 步骤 3: DML 检验（差异甲基化位点）
dmlTest <- DMLtest(BSobj,
  group1 = c("ctrl_rep1", "ctrl_rep2"),
  group2 = c("trt_rep1", "trt_rep2"),
  smoothing = TRUE)

# 查看结果
head(dmlTest)
```

**DMLtest 参数**:
| 参数 | 说明 | 默认值 |
|------|------|--------|
| `smoothing` | 是否使用平滑化 | TRUE |
| `smoothing.span` | 平滑窗口大小 | 500 (bp) |
| `equal.disp` | 假设两组方差相同 | FALSE |
| `BPPARAM` | 并行参数 | SerialParam() |

```r
# 步骤 4: 提取显著 DML
dmls <- callDML(dmlTest,
  p.threshold = 1e-5,
  delta = 0.1)   # 至少 10% 甲基化差异

head(dmls)
cat("Number of significant DMLs:", nrow(dmls), "\n")
```

```r
# 步骤 5: DMR 检测
dmrs <- callDMR(dmlTest,
  p.threshold = 1e-5,
  minlen = 50,        # 最小 DMR 长度 (bp)
  minCG = 3,          # 最少 CG 位点数
  dis.merge = 100,    # 相邻 DML 合并距离 (bp)
  pct.sig = 0.5)      # DML 在 DMR 中的比例

head(dmrs)
cat("Number of DMRs:", nrow(dmrs), "\n")

# 导出 DMR 结果
write.csv(dmrs, "dmr_results_CG.csv", row.names = FALSE)
```

**DMR 参数详解**:
| 参数 | 说明 | 植物建议值 |
|------|------|-----------|
| `p.threshold` | DML p 值阈值 | 1e-5 (可调整) |
| `delta` | 最小甲基化差异 | 0.1 (10%) |
| `minlen` | 最小 DMR 长度 | 50 (bp) |
| `minCG` | 最少 CG 位点 | 3 |
| `dis.merge` | 相邻 DML 合并距离 | 100 (bp) |
| `pct.sig` | DML 最低比例 | 0.5 |

---

## 三种甲基化上下文的 DMR

植物 DNA 甲基化需在三种序列上下文中分别分析：

```r
# CG context DMRs
dmrs_CG <- callDMR(dmlTest_CG, p.threshold = 1e-5)

# CHG context DMRs
dmrs_CHG <- callDMR(dmlTest_CHG, p.threshold = 1e-5)

# CHH context DMRs
dmrs_CHH <- callDMR(dmlTest_CHH, p.threshold = 1e-5)
```

### 上下文特异性 DMR 解读

| 上下文 | DMR 正甲基化 | DMR 负甲基化 | 生物学意义 |
|--------|------------|------------|-----------|
| CG | 超甲基化 | 去甲基化 | 基因沉默/激活 (MET1/CMT3) |
| CHG | 超甲基化 | 去甲基化 | TE 沉默调控 (CMT3) |
| CHH | 超甲基化 | 去甲基化 | RdDM 活性变化 (DRM2) |

---

### 2. methylKit (备选)

**描述**: methylKit 提供更灵活的差显分析和可视化，支持 logistic 回归和 Fisher 精确检验。

```r
library(methylKit)

# 读取数据
file.list <- list("ctrl_rep1_CG.txt", "ctrl_rep2_CG.txt",
                  "trt_rep1_CG.txt", "trt_rep2_CG.txt")

myobj <- methRead(file.list,
  sample.id = list("ctrl1", "ctrl2", "trt1", "trt2"),
  assembly = "genome",
  treatment = c(0, 0, 1, 1),
  context = "CpG")

# 合并样本
meth <- unite(myobj, destrand = FALSE)

# 差异甲基化（Fisher 精确检验）
myDiff <- calculateDiffMeth(meth,
  over.dispersion = "MN",
  adjust = "BH",
  test = "F")

# 筛选 DMR
myDiff_sig <- getMethylDiff(myDiff,
  difference = 25,   # 25% 甲基化差异
  qvalue = 0.01)

write.csv(myDiff_sig, "methylKit_DMR_results.csv")
```

---

## DMR 可视化

### 单个 DMR 信号图
```r
showOneDMR(dmrs[1, ], BSobj,
  ext = 500,   # 延伸 500bp
  main = "DMR Region")
```

### DMR 热图
```r
library(pheatmap)

# 提取 DMR 矩阵
dmr_matrix <- extractMethMatrix(dmrs, BSobj)

# 绘图
pheatmap(dmr_matrix,
  annotation_col = sample_annotation,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  show_rownames = FALSE,
  main = "DMR Methylation Heatmap")
```

---

## 植物 DMR 分析注意事项

1. **区分三种上下文**: 植物必须分别分析 CG、CHG、CHH 的 DMR
2. **重复数要求**: 没有生物学重复时 DSS 仍可运行（使用 `equal.disp = TRUE`），但结果的假阳性率可能升高
3. **测序深度**: CG 甲基化水平比较可靠，CHH 甲基化由于水平低，需要更深的测序深度（>30x）
4. **基因组大小**: 对大型基因组（玉米、小麦），使用染色体/contig分割策略以控制内存
5. **区域类型**: DMR 可能位于基因启动子（调控转录）、基因体（调控剪接）、TE（转座子沉默）、或基因间区（非编码调控）
