# Editing Efficiency Prediction (DeepSpCas9)

**Goal:** Predict the on-target editing efficiency of candidate sgRNAs
**Best for:** Ranking sgRNAs by expected indel frequency and selecting the most active guides

## Prerequisites

- DeepSpCas9 (https://github.com/MyungjaeSong/Paired-Library)
- Python 3.6+ with TensorFlow or PyTorch
- Candidate sgRNA sequences (30bp: 4bp context + 20bp protospacer + 6bp PAM+context)

## DeepSpCas9 Usage

```bash
# Clone and install
git clone https://github.com/MyungjaeSong/Paired-Library.git
cd Paired-Library
pip install -r requirements.txt

# Predict efficiency
python DeepSpCas9.py \
  --input sgRNA_sequences.txt \
  --output efficiency_predictions.csv

# Input format (30bp per line):
# 4bp_upstream + 20bp_protospacer + 3bp_PAM + 3bp_downstream
# Example: ACGTNNNNNNNNNNNNNNNNNNNNNGGNNNN
```

## Python Script for Batch Prediction

```python
import numpy as np
import pandas as pd

def encode_sequence(seq):
    """One-hot encode 30bp sequence for DeepSpCas9"""
    mapping = {'A': [1,0,0,0], 'C': [0,1,0,0],
               'G': [0,0,1,0], 'T': [0,0,0,1]}
    encoded = []
    for nt in seq.upper():
        encoded.extend(mapping.get(nt, [0,0,0,0]))
    return np.array(encoded).reshape(1, 30, 4)

def predict_efficiency(sgRNAs_30bp):
    """Predict editing efficiency for list of sgRNAs"""
    # This requires the trained DeepSpCas9 model
    # Model weights: https://github.com/MyungjaeSong/Paired-Library

    from tensorflow.keras.models import load_model
    model = load_model("DeepSpCas9_model.h5")

    predictions = []
    for seq in sgRNAs_30bp:
        encoded = encode_sequence(seq)
        pred = model.predict(encoded, verbose=0)[0][0]
        predictions.append(pred)

    return predictions

# Usage
candidates_30bp = [
    "ACGT" + sgRNA_20bp + "NGG" + "NNN"
    for sgRNA_20bp in candidate_sgrnas
]
efficiencies = predict_efficiency(candidates_30bp)

# Combine with off-target scores for final ranking
sgRNA_scores = pd.DataFrame({
    "sgRNA": candidate_sgrnas,
    "target_seq": candidates_30bp,
    "deepspcas9_efficiency": efficiencies,
})
sgRNA_scores = sgRNA_scores.sort_values("deepspcas9_efficiency",
                                         ascending=False)
```

## Azimuth (Plant-Specific Models)

```r
# Azimuth 2.0 for plant-specific efficiency prediction
library(azimuth)

# Available organism-specific models:
# - rice_dongjin (Oryza sativa, Dongjin cultivar)
# - arabidopsis (Arabidopsis thaliana, Col-0)

# Predict efficiency for rice
rice_predictions <- predict_efficiency(
  sgRNA_sequences = rice_sgRNAs,
  organism = "rice_dongjin"
)
```

## inDelphi (Repair Outcome Prediction)

```python
# Predict microhomology-mediated repair outcomes
# Web version: https://indelphi.giffordlab.mit.edu/

# For local use: predict frameshift probability
def predict_frameshift_rate(target_seq):
    """
    Predict probability of frameshift-inducing indel
    Based on microhomology score at cut site
    """
    # Cut site is 3bp upstream of PAM
    cut_site = len(target_seq) - 6  # Position of DSB

    # Check for microhomology around cut site
    left_flank = target_seq[max(0, cut_site-10):cut_site]
    right_flank = target_seq[cut_site:min(len(target_seq), cut_site+10)]

    # Microhomology score: number of matching bases after DSB
    mh_score = 0
    for offset in range(1, min(len(left_flank), len(right_flank))):
        if left_flank[-offset:] == right_flank[:offset]:
            mh_score = offset

    # Higher MH score = more predictable deletion size
    # Higher deletion >1bp = better chance of frameshift
    frameshift_prob = min(mh_score / 10.0, 0.8)

    return frameshift_prob
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| DeepSpCas9 score | > 0.5 | At least 50% predicted efficiency |
| Azimuth score | > 0.4 | Plant-specific efficiency threshold |
| Frameshift probability | > 0.6 | High chance of generating frameshift |
| PAM proximity score | Included | Position 18-20 most critical for efficiency |

## Plant-Specific Notes

- Plant chromatin accessibility affects editing efficiency — nucleosome-occupied regions are edited less efficiently
- The efficiency models (DeepSpCas9, Azimuth) were trained primarily on mammalian cells — expect some deviation in plants
- Rice protoplast assays have shown good correlation with DeepSpCas9 predictions (r ~0.6-0.7)
- In planta efficiency is often lower than predicted due to chromatin and delivery limitations
- For plant base editors: different efficiency determinants (use BE-Hive or DeepCBE for CBE, DeepABE for ABE)

## Multi-Factor Efficiency Model

```python
def composite_efficiency_score(sgrna_info):
    """
    Composite score considering multiple factors:
    - DeepSpCas9 prediction
    - GC content
    - Chromatin accessibility (if available)
    - Position in gene (early CDS preferred for KO)
    """
    score = 0.0

    # DeepSpCas9 prediction (weight: 0.5)
    score += 0.5 * sgrna_info["deepspcas9_score"]

    # GC content optimization (weight: 0.2)
    gc = sgrna_info["gc_content"]
    gc_score = max(0, 1 - abs(gc - 50) / 30)  # Optimal around 50%
    score += 0.2 * gc_score

    # Chromatin accessibility (if ATAC-seq/DNase-seq data) (weight: 0.2)
    if "chromatin_score" in sgrna_info:
        score += 0.2 * sgrna_info["chromatin_score"]

    # Position in CDS (weight: 0.1)
    if sgrna_info["rel_position"] < 0.5:  # First half of CDS
        score += 0.1

    return score
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| DeepSpCas9 prediction much higher than in planta | Chromatin, delivery, or species differences | Multiply prediction by 0.4-0.6 for conservative plant estimate |
| Azimuth species mismatch | Using wrong organism model | Switch to correct model or use general DeepSpCas9 |
| inDelphi only works for SpCas9 | Different Cas variants have different repair patterns | Use Cas-specific prediction (Cas12a = different cut pattern) |
| All sgRNAs predicted similar efficiency | Target region has uniform sequence features | Try different target regions (5' UTR, exon 1, exon 2) |
