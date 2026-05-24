# miRNA 靶基因预测 — 工具目录

## 概述

植物 miRNA 靶基因预测是理解 miRNA 功能的关键步骤。植物 miRNA 与靶基因通常具有高度互补性，这一特性被用于靶基因预测算法。

## Key Parameter Decisions

| Parameter | Standard value | When to change | Why |
|-----------|:---:|------|------|
| expectation (max score) | 3.0 | Polyploids: raise to 5.0; highly conserved miRNA families (miR156, miR172): tighten to 2.0 | Polyploid genomes have more potential off-targets from homeologous gene copies; strict cutoff reduces false positives but may miss real targets with tolerated mismatches |
| UPE (unpaired energy) | 25.0 | Non-model species with draft transcriptomes: relax to 30.0 | Incomplete transcriptome assemblies may truncate UTRs, reducing the available pairing context for accurate energy calculation |
| HSP size | 19 | Very short miRNAs (19-20nt): reduce to 17; long miRNAs (23-24nt): increase to 21 | Short miRNAs have fewer contiguous pairing bases; HSP length must match the seed pairing window of the query miRNA |
| translation inhibition range | 9-11 nt | Non-canonical miRNA binding with central mismatch: adjust to 10-13 nt | Plant miRNAs can repress translation rather than cleave when central mismatches are present; the cleavage-competent pairing window shifts accordingly |

## 推荐工具

### 1. psRNATarget (强烈推荐)

**描述**: 专门为植物设计的 miRNA 靶基因预测工具，基于 Smith-Waterman 算法和改进的评分方案。

**关键特性**:
- 基于植物 miRNA-mRNA 高度互补的特性
- 支持大规模批量预测
- 提供靶位点翻译抑制评估

**输入**:
- miRNA 序列 (FASTA)
- 转录本/转录组序列 (FASTA)
- 用户自定义罚分参数

**输出**:
- 靶基因列表及互补性评分
- 靶位点详细比对信息

**参考**: http://plantgrn.noble.org/psRNATarget/

### 2. TargetFinder (备选)

**描述**: 另一种植物 miRNA 靶基因预测工具，适用于棉花等特定物种。

**特性**:
- 基于降解组数据验证的一致性
- 支持自定义评分矩阵

**平台**: 网页服务 / 本地运行

## Plant-Specific Notes

- **High complementarity requirement**: Unlike animal miRNAs, plant miRNAs typically require near-perfect complementarity (<=4 mismatches) for target cleavage. Use psRNATarget with stricter settings than animal target prediction tools like miRanda or TargetScan — the expectation parameter should be <=5.0 for plant data.
- **Degradome-seq validation**: psRNATarget predictions should be validated with degradome-seq (PARE) data when available. Cleaved miRNA targets show a characteristic peak at the 10th-11th nucleotide of the miRNA binding site in degradome signal. Public plant degradome datasets are available for Arabidopsis, rice, maize, soybean, and tomato.
- **Polyploid target space complexity**: In polyploid species, a single miRNA may target homeologous gene copies with slightly different complementarity due to sequence divergence between subgenomes. Analyze each subgenome separately when possible, and compare target gene ontology enrichment between homeologous target sets to identify conserved vs diverged regulatory functions.
- **miRNA-target regulatory modules**: Plant miRNAs often regulate multiple genes in the same pathway (e.g., miR156-SPL phase transition, miR319-TCP leaf development, miR164-NAC organ boundary). When a miRNA targets multiple members of a gene family, validate the regulatory network topology with co-expression analysis or dual-luciferase assays.
- **Non-canonical targets**: Some plant miRNAs bind with central mismatches and repress translation rather than trigger cleavage. psRNATarget's translation inhibition score helps identify these cases, but they require experimental validation via proteomics or reporter assays.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| No targets found for known miRNA | Expectation too strict or transcriptome missing UTRs | Relax expectation to 5.0; verify transcriptome assembly includes full-length cDNAs with UTRs |
| Thousands of targets per miRNA | Expectation too loose or short transcriptome with many false matches | Tighten expectation to 2.0; filter by target accessibility (consider secondary structure at binding site) |
| psRNATarget server timeout or XML error | Poor internet connection or server overload during batch submission | Use local psRNATarget installation; submit jobs in groups of 50 sequences or fewer |
| False positive targets inflated in polyploids | Homeologous transcripts treated as independent targets inflates target count | Collapse homeologous transcripts (CD-HIT at 98% identity) before prediction, or flag redundant targets post hoc |
