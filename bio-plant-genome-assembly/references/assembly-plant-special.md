# -- 

## 1. 

### 1.1 




|  |  |  | hifiasm  | Flye  |
|------|-----------|------|--------------|---------|
|  | ~135 Mb |  | 8 GB | 4 GB |
|  | ~430 Mb |  | 16 GB | 8 GB |
|  | ~2.4 Gb |  | 64 GB | 32 GB |
|  | ~1.1 Gb |  | 32 GB | 16 GB |
|  | ~17 Gb |  | 256 GB+ | 128 GB+ |

### 1.2 

GenomeScope 2.0k-mer:

- **k-mer**k=21 
- ****
- ****
- ****0.5%; >1% hifiasm  Purge_dups

### 1.3 



|  | HiFi  | ONT  | Illumina  |
|------|--------|---------|-----------|
|  (<500Mb) | 20-30x | 30-50x | 30-50x |
|  (500Mb-3Gb) | 25-40x | 40-60x | 40-60x |
|  (>3Gb) | 30-50x | 50-80x | 50-80x |

---

## 2. 

### 2.1 

|  |  | 
|------|------|---------|
| HiFi (PacBio Revio/Sequel IIe) | > 99.9% | hifiasm  |
| ONT (R10.4.1) | 98-99% (Q20+) | Flye + Medaka |
| ONT (R9.4.1) | 90-95% | Flye + Medaka + Pilon |
| + ONT  |  | hifiasm `--ul` |

### 2.2 

hifiasm HiFi; Flye ONT

### 2.3 hifiasm 

- ****`*.bp.p_ctg.gfa`
- ****hifiasm `*.bp.hap1.p_ctg.gfa` `*.bp.hap2.p_ctg.gfa`
- **purge_dups**> 1% 

---

## 3. 

### 3.1 



|  |  |
|------|------|
| ONT  | Medaka |
| PacBio HiFi  | gcpp |
| +  | Pilon |
|  | Merqury QV > 30 |

### 3.2 

- **QV**-10 x log10() 30=99.9%; 40=99.99%
- **QV > 40**; QV > 30

---

## 4. 

### 4.1 

|  |  |
|---------|---------|
| embrophyt a_odb10 |  |
| eudicots_odb10 |  |
| liliopsida_odb10 |  |
| brassicales_odb10 |  |
| poales_odb10 |  |
| fabales_odb10 |  |

### 4.2 

|  |  |  | 
|------|------|------|------|
| Complete © | > 95% | > 90% | > 80% |
| Single (S) | > 70% | > 60% | > 50% |
| Duplicated (D) | < 5% | < 10% | < 15% |

---

## 5. Hi-C 

### 5.1 

|  |  | 
|------|------|------|
| YAHS | Contig +  | 
| SALSA2 | Contig +  | 
| 3D-DNA/Juicer |  | 

### 5.2 Hi-C 

- DpnII (GATC)
- Hi-C 
- Hi-C (100-500kb)

---

## 6. 

### 6.1 

1. 
2. 
3. 
4.  (TE) 
5.  rDNA 

### 6.2 

- **HiFi **hifiasm `*.bp.hap1`  `*.bp.hap2`
- **ONT **Flye Medaka
- **+ **Pilon 

### 6.3 

- Hi-C Hi-C contig 
- 

---

## 7. 

### 7.1 

|  |  |
|------|------|
| HiFi  | hifiasm |
| ONT  | Flye + Medaka |
| +  | MaSuRCA   |
| Hi-C  | YAHS  |
|  | BUSCO + Merqury |

### 7.2 

4 :

1.  +  + 
2. hifiasm  Flye (+ purge_dups )
3. Medaka/gcpp/Pilon 
4. BUSCO + Merqury  + Hi-C contact map

---

## 8. 

### hifiasm "Killed"

- hifiasm > 1 Gb  64 GB RAM
- `-s`  ( )
- ; 

### Flye 

- ; > 20x 
- `--genome-size` 

### BUSCO Duplicated > 20%

- purge_dups 
- 
- WGD (Whole Genome Duplication)

### Hi-C 

- Hi-C contact map 15%
-  (juicebox)
- Hi-C 300M+ valid pairs

---

## 9. 

- hifiasm: https://github.com/chhylp123/hifiasm
- Flye: https://github.com/fenderglass/Flye
- BUSCO: https://busco.ezlab.org/
- Merqury: https://github.com/marbl/merqury
- YAHS: https://github.com/c-zhou/yahs
- GenomeScope 2.0: http://qb.cshl.edu/genomescope/genomescope2.0/
- purge_dups: https://github.com/dfguan/purge_dups
- Plant genome assembly review: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02587-y
