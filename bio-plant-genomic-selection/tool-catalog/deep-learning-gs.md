# Deep Learning for Genomic Selection

**Goal:** Capture non-additive and epistatic effects using neural networks
**Best for:** Very large populations (≥2000) with high-density markers (≥100K)

## Prerequisites
- Python 3.10+, PyTorch or TensorFlow, scikit-learn
- Large training population
- GPU recommended

## Basic MLP Approach

```python
import torch
import torch.nn as nn
from sklearn.model_selection import KFold

class GS_MLP(nn.Module):
    def __init__(self, n_snps, hidden=[256, 128, 64]):
        super().__init__()
        layers = []
        prev = n_snps
        for h in hidden:
            layers.extend([nn.Linear(prev, h), nn.ReLU(), nn.Dropout(0.2)])
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)
```

## When DL Might Help

- >2000 training samples with >100K markers
- Traits with known epistatic/interaction effects
- Multi-omics integration (genotype + transcriptome + metabolome)

## When DL Probably Won't Help

- <1000 samples — overfitting risk is too high
- Highly inbred crops — additive effects dominate, GBLUP is near optimal
- Single environment, moderate heritability — non-additive variance is small

## Plant-Specific Notes

- DL for GS in plants is still experimental — GBLUP/Bayes remain the production methods
- For hybrid prediction (combining ability), DL may capture dominance effects that GBLUP misses
- Always benchmark DL against GBLUP — if DL isn't >5% better, use GBLUP
