# ATAC-seq -- 

## 1. 

### 1.1 

- **true**

`inspect_data.sh`  `data_profile.json`



 `species_detected`  `breeding_system` 
 `bio-plant-infra/references/species-cheatsheet.md` 


chr1..chr10  Chr01..Chr12  1A..7D 
ATAC-seq 

### 1.2 

- **peak 30,000-100,000**ATAC-seq 
- **peak 10,000-30,000**
- **peak < 10,000**

### 1.3 ATAC-seq 

Tn5 ATAC-seq 

- ****Tn5 
- ****
- ****

### 1.4 

- **n > 2**DiffBind + DESeq2
- **n = 2**DESeq2 (FDR )
- **n = 1**

### 1.5 

- ****NFR (Nucleosome-Free Region)
- ****Mono-Nucleosome
- ****Di-Nucleosome

Tn5 TSS 

---

## 2. Peak Calling

### 2.1 MACS2 ATAC-seq 

MACS2 `--nomodel --shift -100 --extsize 200`Tn5 

### 2.2 

- **--nomodel******
- **-q 0.05**FDR
- ****ATAC-seq `--keep-dup all`

### 2.3 

- 5-20M reads 
- 20-40M reads 
- 40-80M reads 

---

## 3. 

### 3.1 

|  |  |  |
|------|------|------|
| Promoter (<=1kb) |  |  |
| Distal intergenic |  |  |
| Intron |  |  |
| Exon |  |  |

### 3.2 

ChIPseeker + plant OrgDb 

---

## 4. 

### 4.1 

1. **+ **2 
2. DESeq2 (DiffBind  DESeq2)
3. 

### 4.2 

- FDR < 0.05
- |log2FC| > 0 (DiffBind )

---

## 5. TF Footprinting

### 5.1 

Tn5 TF 

### 5.2 TOBIAS 

1. ATACorrect: Tn5 
2. ScoreMotifs: motif footprint 
3. BINDetect: TF 

### 5.3 

- ATAC-seq  50M reads
- TF 

---

## 6. Motif 

### 6.1 

|  | URL |
|------|-----|
| PlantTFDB | http://planttfdb.gao-lab.org/ |
| JASPAR Plants | https://jaspar.genereg.net/ |
| CIS-BP | http://cisbp.ccbr.utoronto.ca/ |
| PlantPAN | http://plantpan.itps.ncku.edu.tw/ |

### 6.2 

- De novo motif: HOMER  MEME-ChIP
- : AME (MEME Suite)
- Tomtom  de novo motif  TF

---

## 7. 

### 7.1 

- ATAC-seq 
- MACS2  `--nomodel --shift -100 --extsize 200`
- 2  DiffBind
-  TOBIAS

### 7.2 

3 :

1. 

2. 
3. 

---

## 8. 

### 

- 50%+ cpDNA 
- PCR  Tn5 
- ATAC-seq ChIP-seq Input 
- Bowtie2 `--very-sensitive`-X 2000

### 

- Tn5 
- 50M reads
-  ( 2x)hifiasm 
- TOBIAS --split

### ChIP-eq ATAC-seq 

- **Input **ChIP-seq Input DNATn5 
- **Peak **ChIP-seq narrowATAC-seq MACS2 `--nomodel --shift -100`
- **FRiP **ChIP TF > 5%ATAC-seq > 30%
- **Footprinting**ATAC-seq ; ChIP-seq 

---

## 9. 

- MACS2: https://github.com/macs3-project/MACS
- deepTools: https://deeptools.readthedocs.io/
- ChIPseeker: https://bioconductor.org/packages/ChIPseeker/
- DiffBind: https://bioconductor.org/packages/DiffBind/
- TOBIAS: https://github.com/loosolab/TOBIAS
- HOMER: http://homer.ucsd.edu/homer/
- MEME Suite: https://meme-suite.org/meme/
- ENCODE ATAC-seq Guidelines: https://www.encodeproject.org/atac-seq/
