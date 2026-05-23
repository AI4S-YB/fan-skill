# 甲基化-表达关联分析 -- 工具目录

## 概述

将 DNA 甲基化数据与基因表达数据（RNA-seq）整合，分析表观遗传调控对转录的影响。植物中需要关注启动子甲基化和基因体甲基化（gbM）对表达的不同效应。

## 推荐分析策略

### 1. 启动子甲基化-表达关联

**描述**: 基因启动子附近（TSS +/- 2kb）的 DNA 甲基化通常与基因表达负相关。

#### 数据准备

```r
# 甲基化数据
library(methylKit)
meth_promoter <- regionCounts(myobj, promoters)
promoter_meth <- as.data.frame(meth_promoter)

# 表达数据 (来自 RNA-seq: FPKM/TPM 矩阵)
expression <- read.csv("gene_expression_tpm.csv",
  row.names = 1)
```

#### 关联分析

```r
# 计算每个基因的启动子甲基化与表达的相关性
cor_results <- data.frame(
  gene = character(),
  cor_coef = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in intersect(rownames(promoter_meth),
                       rownames(expression))) {
  meth_values <- as.numeric(promoter_meth[gene, ])
  expr_values <- as.numeric(expression[gene, ])

  if (sd(meth_values) > 0 && sd(expr_values) > 0) {
    cor_test <- cor.test(meth_values, expr_values,
      method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      gene = gene,
      cor_coef = cor_test$estimate,
      p_value = cor_test$p.value
    ))
  }
}

# 多重检验校正
cor_results$padj <- p.adjust(cor_results$p_value,
  method = "BH")

# 显著负相关的基因 (启动子高甲基化 -> 低表达)
neg_cor_genes <- cor_results[
  cor_results$cor_coef < -0.5 & cor_results$padj < 0.05,
]
cat("Genes with significant negative correlation:",
  nrow(neg_cor_genes), "\n")
```

### 2. 甲基化与表达的分组比较

**在有无 DMR 的条件下比较表达水平**:

```r
# 将基因分为两组：有 DMR 在启动子 vs 无 DMR
dmr_genes <- unique(dmrs$gene)
non_dmr_genes <- setdiff(all_genes, dmr_genes)

# 比较两组表达
dmr_expr <- expression[dmr_genes, ]
non_dmr_expr <- expression[non_dmr_genes, ]

# Wilcoxon 检验
wilcox_test <- wilcox.test(
  rowMeans(dmr_expr),
  rowMeans(non_dmr_expr)
)

# 可视化
boxplot(list(
  "With DMR" = rowMeans(dmr_expr),
  "Without DMR" = rowMeans(non_dmr_expr)
), main = "Expression vs Promoter DMR",
col = c("salmon", "steelblue"),
ylab = "Mean Expression (TPM)")
```

### 3. 基因体甲基化 (gbM) 分析

```r
# 将基因分为 gbM 阳性/阴性
gbM_positive <- gbM_genes  # CG 甲基化基因体 > 阈值
gbM_negative <- setdiff(all_genes, gbM_positive)

# 比较两组的表达水平分布
gbm_expr_pos <- expression[gbM_positive, ]
gbm_expr_neg <- expression[gbM_negative, ]

# gbM 基因通常有中等稳定的表达
boxplot(list(
  "gbM+" = rowMeans(gbm_expr_pos),
  "gbM-" = rowMeans(gbm_expr_neg)
), main = "Expression by gbM Status",
col = c("orange", "gray"),
ylab = "Mean Expression (TPM)")
```

### 4. 多组学整合可视化

```r
library(circlize)
library(ComplexHeatmap)

# 联合热图：甲基化 + 表达
# 构建联合矩阵
combined_matrix <- cbind(
  scale(t(methylation_sig_dmrs)),
  scale(t(expression_sig_genes))
)

# 分层列注释
ha_col <- HeatmapAnnotation(
  DataType = c(rep("Methylation", n_sig_dmrs),
               rep("Expression", n_sig_genes)),
  col = list(DataType = c("Methylation" = "blue",
                          "Expression" = "red"))
)

# 联合热图
Heatmap(combined_matrix,
  top_annotation = ha_col,
  name = "Z-score",
  show_row_names = FALSE,
  cluster_columns = FALSE,
  column_split = ha_col$DataType)
```

### 5. 四象限分析图

```r
# 四象限图: 甲基化变化 vs 表达变化
# x 轴: DMR 甲基化差异 (delta beta)
# y 轴: 对应基因表达 log2FC

plot(dmr_delta[, "CG"],
  expr_log2FC[common_genes],
  xlab = "Methylation Difference (Delta Beta)",
  ylab = "Expression log2 Fold Change",
  main = "Methylation-Expression Relationship",
  pch = 16, col = rgb(0, 0, 0, 0.3))

abline(h = 0, v = 0, lty = 2, col = "gray")

# 四象限统计
Q1 <- sum(dmr_delta > 0 & expr_log2FC > 0)  # 高甲基化+高表达
Q2 <- sum(dmr_delta < 0 & expr_log2FC > 0)  # 低甲基化+高表达
Q3 <- sum(dmr_delta < 0 & expr_log2FC < 0)  # 低甲基化+低表达
Q4 <- sum(dmr_delta > 0 & expr_log2FC < 0)  # 高甲基化+低表达

cat("Q1 (Hyper-Up):", Q1, "\n")
cat("Q2 (Hypo-Up):", Q2, "\n")
cat("Q3 (Hypo-Down):", Q3, "\n")
cat("Q4 (Hyper-Down):", Q4, "\n")
```

---

## 植物特有甲基化-表达关系

### 启动子甲基化的负调控

在植物中，启动子甲基化通常抑制转录：
- **CG 启动子甲基化**: 与转录沉默强相关
- **CHH 启动子甲基化**: RdDM 介导的沉默（可遗传）
- **CHG 启动子甲基化**: CMT3 维持的异染色质标记

### 基因体甲基化 (gbM) 的正/中性效应

**基因体甲基化** (gene body methylation, gbM) 是植物的独特特征：
- 存在于约 20-30% 的拟南芥基因
- 是 CG 特异性的（非 CHG 或 CHH）
- 排除启动子和转录终止位点
- gbM 基因通常是组成型表达的管家基因
- gbM 基因的表达水平比非 gbM 基因更稳定

### 转座子甲基化的表达效应

植物中 TE 甲基化通过两种机制影响邻接基因的表达：
1. **TE 插入在启动子**: 招募 DNA 甲基化，抑制邻接基因
2. **TE 插入在 intron**: 可能影响剪接（splicing）

---

## 常用统计方法

| 方法 | 适用场景 |
|------|---------|
| Spearman 相关 | 有序但非线性的甲基化-表达关系 |
| Wilcoxon 秩和检验 | 两组比较（有/无 DMR） |
| 线性回归 | 假设线性甲基化-表达关系 |
| 广义加性模型 (GAM) | 非线性甲基化-表达关系 |
