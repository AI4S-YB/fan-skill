# miRNA 预测 — 工具目录

## 概述

miRNA 预测是植物小 RNA 分析的核心步骤。本目录涵盖基于测序数据的已知 miRNA 鉴定和新 miRNA 发现。

## 推荐工具

### 1. miRDeep2 (推荐)

**描述**: 基于二级结构和测序 reads 的 miRNA 预测工具。

**适用场景**:
- 模式物种：结合 miRBase 已知 miRNA 进行鉴定
- 非模式物种：从头预测新 miRNA

**关键参数**:
- `-g <genome>`: 参考基因组 FASTA
- `-f <fasta>`: 测序 reads (FASTA 格式)
- `-d <mirbase>`: miRBase 已知 miRNA (模式物种)
- `-t <species>`: 物种名称

**输出**:
- `result.csv`: 预测的 miRNA 列表及其评分
- `mirna_structure.pdf`: miRNA 前体二级结构

**参考**: https://github.com/rajewsky-lab/mirdeep2

### 2. ShortStack (备选)

**描述**: 植物小 RNA 分析的综合工具，支持 miRNA、siRNA、phasiRNA 鉴定。

**适用场景**: 需要同时分析多种小 RNA 类型时使用。

**输出**:
- miRNA 注释
- siRNA/phasiRNA 簇
- 小 RNA 比对统计
