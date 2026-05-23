# -- 

## 1. 

### 1.1 

- **CG (CpG)**MET1 
- **CHG**CMT3 
- **CHH**DRM2 (RdDM )


|  | CG (%) | CHG (%) | CHH (%) |
|------|--------|---------|---------|
|  | 24-30 | 6-10 | 2-5 |
|  | 30-45 | 10-15 | 3-8 |
|  | 50-60 | 20-30 | 2-5 |
|  | 40-55 | 15-25 | 3-6 |
|  | 50-65 | 20-30 | 5-10 |

### 1.2 

- **WGBS**
- **RRBS**MspI 
- **EM-seq**

### 1.3 

`inspect_data.sh`bisulfite conversion rate

---

## 2. Bismark 

### 2.1 

**WGBS `--non_directional`**Bismark 

Bisulfite  C>T  G>A 

### 2.2 WGBS vs RRBS 

|  | WGBS | RRBS |
|------|------|------|
|  | `--non_directional` | `--directional` |
|  | 15-30x | 10-20x |
|  |  | CpG  |
| Bismark  |  | `--rrbs` |

### 2.3 

- Bismark 40-60%bisulfite 
- 4  `--non_directional`

---

## 3. 

### 3.1 CG/CHG/CHH 

DSS  methylKit

### 3.2 

- **CG **MET1 
- **CHG **RdDM 
- **CHH **DRM2/RdDM 

### 3.3 DMR 

DSS  DMRp.threshold = 1e-5

|  | DMR  | DMR  |
|------|------|------|
| CG |  () |  |
| CHG |  |  |
| CHH | RdDM  | RdDM  |

---

## 4. 

### 4.1 

- WGBS 30xCHH 
- RRBS CpG CHG/CHH 
- 5x 

### 4.2 

- **CpG**~200bp 
- **CHG**: CpG 
- **CHH**: 5-10% 

---

## 5. 

### 5.1 

- **gbM (gene body methylation)**CG 20-30%
- gbM 

### 5.2 TE 

- TE  CG/CHG/CHH 
- RdDM  24-nt siRNA  TE  CHH 
- TE 

### 5.3 

- **RdDM (RNA-directed DNA Methylation)** 24nt siRNA 
- **MET1/CMT3/DRM2** 
-  (imprinting)  (endosperm demethylation)

---

## 6. 




---

## 7. 

### 7.1 

1. WGBS  Bismark `--non_directional`
2. Bismark_met_extract  CX_context
3. CG/CHG/CHH 
4. DSS  2 
5. 

### 7.2 

3 :

1.  CG/CHG/CHH  WGBS 
2. 
3. 

---

## 8. 

### WGBS 

- 30x CHH 
- WGBS `--non_directional`
- CG  10x  coverage2cytosine

### 

-  (T->C)
- CpG  (CpG > 10x)
- WGBS  RRBS  DMR

### 

-  CG  (gbM)
- WGBS CHH  (RdDM)
- WGBS  C>T 

---

## 9. 

- Plant MethylDB: http://epigenome.genetics.uga.edu/PlantMethylDB/
- MethBank: https://ngdc.cncb.ac.cn/methbank/
- Plant Epigenome Database: http://systemsbiology.cau.edu.cn/chromstates/
- RdDM pathways: https://www.nature.com/articles/nrg3885
- plant gbM review: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-x
