# CRISPR Design Visualization

**Goal:** Create visualizations of sgRNA target sites, vector maps, and editing outcomes
**Best for:** Communicating CRISPR design to collaborators and for publication figures

## sgRNA Target Map

```python
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_sgrna_target_map(gene_structure, sgRNAs, target_gene_name):
    """
    Plot gene structure with sgRNA target positions

    gene_structure: dict with exons, CDS, UTR coordinates
    sgRNAs: list of dicts with position, strand, sequence, scores
    """
    fig, ax = plt.subplots(figsize=(12, 4))

    gene_length = gene_structure['end'] - gene_structure['start']

    # Draw gene body line
    ax.plot([0, gene_length], [0, 0], 'k-', linewidth=2)

    # Draw exons (boxes)
    for exon in gene_structure['exons']:
        exon_start = exon['start'] - gene_structure['start']
        exon_end = exon['end'] - gene_structure['start']
        if exon['type'] == 'CDS':
            color = '#2c3e50'
            height = 0.3
        else:  # UTR
            color = '#95a5a6'
            height = 0.15
        rect = patches.Rectangle((exon_start, -height/2),
                                 exon_end - exon_start, height,
                                 linewidth=1, edgecolor='black',
                                 facecolor=color)
        ax.add_patch(rect)

    # Draw sgRNA targets (arrows)
    for i, sgrna in enumerate(sgRNAs):
        pos = sgrna['position'] - gene_structure['start']
        y_offset = 0.4 + i * 0.15
        ax.annotate('', xy=(pos, y_offset),
                    xytext=(pos - 20, y_offset),
                    arrowprops=dict(arrowstyle='->', color='#e74c3c', lw=2))
        ax.text(pos + 10, y_offset,
                f"{sgrna['id']} (Eff:{sgrna['efficiency']:.0f})",
                fontsize=8, va='center')

    ax.set_xlim(0, gene_length)
    ax.set_ylim(-0.5, 0.5 + len(sgRNAs) * 0.15)
    ax.set_xlabel("Position (bp)")
    ax.set_yticks([])
    ax.set_title(f"sgRNA Target Map — {target_gene_name}")

    plt.tight_layout()
    plt.savefig("sgrna_target_map.png", dpi=300, bbox_inches='tight')
```

## Off-Target Circos Plot

```python
def plot_off_target_circos(on_target_chr, on_target_pos, off_targets_df):
    """
    Circos-style plot showing on-target and off-target positions
    across the genome

    Requires: chromosome lengths file
    """
    from circlize import Circos

    circos = Circos(chromosome_sizes)

    # Draw chromosomes
    for chr_name, chr_len in chromosome_sizes.items():
        circos.add_sector(chr_name, chr_len)

    # Mark on-target site (green)
    circos.add_point(on_target_chr, on_target_pos,
                     color='#2ecc71', size=8, label='On-target')

    # Mark off-target sites (colored by mismatch count)
    for _, row in off_targets_df.iterrows():
        if row['Mismatch'] == 0:
            color = '#e74c3c'  # Red: 0 mismatches (high risk)
        elif row['Mismatch'] == 1:
            color = '#f39c12'  # Orange: 1 mismatch
        else:
            color = '#3498db'  # Blue: 2+ mismatches

        circos.add_point(row['Chr'], row['Pos'],
                         color=color, size=5)

    circos.draw()
    plt.savefig("off_target_circos.png", dpi=300)
```

## Editing Efficiency Bar Plot

```r
library(ggplot2)

# Efficiency comparison of candidate sgRNAs
sgrna_scores <- data.frame(
  sgRNA = c("sgRNA1", "sgRNA2", "sgRNA3", "sgRNA4", "sgRNA5"),
  Efficiency = c(0.72, 0.65, 0.58, 0.45, 0.38),
  Specificity = c(0.85, 0.92, 0.78, 0.95, 0.88)
)

ggplot(sgrna_scores, aes(x = reorder(sgRNA, -Efficiency))) +
  geom_col(aes(y = Efficiency), fill = "#3498db", alpha = 0.8) +
  geom_point(aes(y = Specificity), color = "#e74c3c", size = 4) +
  geom_line(aes(y = Specificity, group = 1), color = "#e74c3c") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     name = "Score",
                     sec.axis = sec_axis(~., name = "Specificity (CFD)")) +
  labs(x = "", title = "sgRNA Efficiency and Specificity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

## HDR Donor Map

```python
def plot_hdr_donor_map(donor_design):
    """
    Visualize HDR donor template structure
    """
    fig, ax = plt.subplots(figsize=(10, 2))

    y = 0
    x = 0

    # Left homology arm
    ax.add_patch(patches.Rectangle((x, -0.2), donor_design['left_arm_len']/100,
                                   0.4, facecolor='#3498db', alpha=0.7))
    ax.text(x + donor_design['left_arm_len']/200, 0, "Left HA",
            ha='center', va='center', fontsize=10, color='white')
    x += donor_design['left_arm_len']/100

    # Insertion
    ax.add_patch(patches.Rectangle((x, -0.2), donor_design['insert_len']/100,
                                   0.4, facecolor='#e74c3c', alpha=0.7))
    ax.text(x + donor_design['insert_len']/200, 0, "Insert",
            ha='center', va='center', fontsize=10, color='white')
    x += donor_design['insert_len']/100

    # Right homology arm
    ax.add_patch(patches.Rectangle((x, -0.2), donor_design['right_arm_len']/100,
                                   0.4, facecolor='#3498db', alpha=0.7))
    ax.text(x + donor_design['right_arm_len']/200, 0, "Right HA",
            ha='center', va='center', fontsize=10, color='white')

    # Mark mutated PAM
    if donor_design.get('pam_pos'):
        ax.axvline(x=donor_design['pam_pos']/100, color='#f39c12',
                   linestyle='--', linewidth=2)
        ax.text(donor_design['pam_pos']/100, 0.3, "PAM mutated",
                ha='center', fontsize=8, color='#f39c12')

    ax.set_xlim(0, x + 1)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel("Position (x100 bp)")
    ax.set_yticks([])
    ax.set_title("HDR Donor Template Design")

    plt.tight_layout()
    plt.savefig("hdr_donor_map.png", dpi=300)
```

## Plant-Specific Figure Notes

- Multi-panel figure for publication: (a) sgRNA target map on gene, (b) off-target Circos, (c) editing efficiency bar chart, (d) vector map
- Show the plant species and Cas variant used in the title
- For multiplex editing: show all target genes in a single figure with connecting arrows
- For base editing: show the editing window on the sequence (positions 4-8 highlighted)
- Include scale bars on all genomic maps

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| Gene structure plot too crowded | Too many sgRNAs or large introns | Collapse introns, use "//" notation |
| Circos plot unreadable | Too many off-targets | Filter to <=2 mismatches or top 50 |
| Sequence too small to read | Font size too small for publication | Export as PDF (vector), use 14pt font |
| Colors not colorblind-friendly | Red/green only | Use viridis or colorblind-friendly palette |
