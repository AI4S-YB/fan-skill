# HDR Template Design

**Goal:** Design donor templates for homology-directed repair (HDR) knock-in experiments in plants
**Best for:** Precise gene insertion, epitope tagging, and allele replacement via CRISPR/Cas9

## Prerequisites

- sgRNA target sequence and cut site coordinates
- Insertion sequence (gene, tag, reporter)
- Target locus sequence (FASTA, including flanking regions)

## HDR Donor Template Design

### Basic Design Rules

```
Left Homology Arm --- Insertion Sequence --- Right Homology Arm
    400-800bp         (gene/tag/marker)         400-800bp

Key features:
1. Cut site positioned at junction of homology arms
2. PAM sequence MUTATED (silent mutation) in donor to prevent recutting
3. sgRNA target sequence MUTATED in donor (prevent donor cleavage)
```

### Donor Template Construction

```python
def design_hdr_donor(target_locus, insertion_seq, cut_site_pos,
                     left_arm_len=600, right_arm_len=600):
    """
    Design HDR donor template for CRISPR knock-in

    Parameters:
    - target_locus: str, genomic DNA sequence around target site
    - insertion_seq: str, sequence to insert
    - cut_site_pos: int, position of Cas9 cut within target_locus
    - left_arm_len: length of left homology arm
    - right_arm_len: length of right homology arm
    """
    # Extract homology arms
    left_arm = target_locus[max(0, cut_site_pos - left_arm_len):cut_site_pos]
    right_arm = target_locus[cut_site_pos:min(len(target_locus), cut_site_pos + right_arm_len)]

    # Mutate PAM in right arm (if PAM is in right arm)
    # NGG -> mutate one G to prevent Cas9 binding
    pam_pos = target_locus.find("GG", cut_site_pos - 3)
    if pam_pos != -1:
        # Introduce silent mutation: G -> A at wobble position
        right_arm_mutated = right_arm[:pam_pos - cut_site_pos] + \
                            right_arm[pam_pos - cut_site_pos:].replace("GG", "GA", 1)

    # Assemble donor
    donor = left_arm + insertion_seq + right_arm_mutated

    # Add restriction sites for cloning (optional)
    # donor = "GAATTC" + donor + "AAGCTT"

    return {
        "donor_sequence": donor,
        "left_arm": left_arm,
        "right_arm": right_arm_mutated,
        "insertion": insertion_seq,
        "total_length": len(donor),
        "pam_mutated": pam_pos != -1
    }
```

## Cloning Strategy

### Golden Gate Assembly

```python
# Design primers for Golden Gate assembly
def design_golden_gate_primers(donor_parts, vector_backbone):
    """
    Design BsaI/BbsI Golden Gate assembly for donor cloning
    """
    primers = []
    overhangs = []

    # Level 0: individual parts with BsaI sites
    for i, part in enumerate(donor_parts):
        # Add BsaI sites + unique 4bp overhang
        left_primer = f"gaggtctcA{overhangs[i]}{part[:20]}"
        right_primer = f"cgtctcacT{overhangs[i+1]}{part[-20:][::-1]}"
        primers.append((left_primer, right_primer))

    return primers, overhangs
```

## Plant-Specific HDR Strategies

### Geminivirus Replicon System

The most effective HDR approach for plants:

```python
# Design geminivirus-based donor replicon
def design_geminivirus_donor(donor_template, rep_gene_source="Wheat_Dwarf_Virus"):
    """
    Design donor template for geminivirus replicon amplification

    Key components:
    - LIR (Long Intergenic Region): origin of replication
    - SIR (Short Intergenic Region): required for replication
    - Rep/RepA: replication initiator proteins
    - Donor template between LIR elements
    """

    replicon = {
        "LIR_left": "...",    # WDV LIR sequence
        "donor_template": donor_template,
        "LIR_right": "...",   # WDV LIR sequence
        "Rep_expression_cassette": "35S:Rep",
    }

    return replicon
```

## PEG-RNA Design (Prime Editing)

```python
def design_pegRNA(target_seq, edit_type, edit_sequence, position):
    """
    Design prime editing guide RNA (pegRNA)

    Components:
    - spacer (20nt sgRNA protospacer)
    - scaffold (sgRNA scaffold)
    - RTT (reverse transcription template): the edit + homology
    - PBS (primer binding site): ~13nt binding to nicked DNA
    """

    # Protospacer + scaffold
    spacer = target_seq[:20]  # sgRNA spacer
    scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

    # PBS design (~13 nt)
    # PBS binds to the 3' end of the nicked DNA strand
    pbs_len = 13
    nick_site_flank = target_seq[20:20+pbs_len]  # Sequence after cut
    pbs = str(Seq(nick_site_flank).reverse_complement())

    # RTT design
    # Contains: edit + homology to template after edit
    rtt_len = 15  # RTT recommended length for small edits
    if edit_type == "insertion":
        rtt = edit_sequence + target_seq[20:20+rtt_len]
    elif edit_type == "substitution":
        # Target position within protospacer or RTT
        rtt = target_seq[position:position+len(edit_sequence)].replace(
            target_seq[position:position+len(edit_sequence)],
            edit_sequence
        ) + target_seq[position+len(edit_sequence):position+len(edit_sequence)+rtt_len]

    # Full pegRNA
    pegRNA = spacer + scaffold + rtt + pbs

    return pegRNA
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| Homology arm length | 400-800 bp | Longer = higher HDR efficiency |
| Total donor length | < 2,000 bp | Cloning and delivery limitations |
| PAM mutation | Silent or intronic | Prevent donor re-cutting |
| HDR efficiency (plant) | 0.5-5% | Low baseline, use enrichment strategies |
| Prime editing efficiency (plant) | 1-10% | Higher than HDR for small edits |

## Plant-Specific Notes

- HDR efficiency in plants is very low (<5%) without enhancement strategies
- Geminivirus replicon can increase HDR by 10-100x (amplifying donor template in planta)
- Cell cycle synchronization to S/G2 phase improves HDR (HDR is active during these phases)
- For small edits (<40 bp), prime editing is now preferred over HDR
- T0 plants are chimeric — screen T1 generation for clean knock-in events
- Marker-free selection: co-deliver editing components on separate T-DNAs that segregate in T2

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| No knock-in events in T0 | HDR efficiency too low | Use geminivirus replicon or prime editing |
| Donor integrated at random locus | Donor integrated via NHEJ | Add NHEJ inhibitor or increase homology arm length |
| Donor recut by Cas9 | PAM not mutated in donor | Mutate PAM sequence (silent mutation) |
| pegRNA design errors | Incorrect PBS or RTT length | Use pegFinder (prediction tool), test 2-3 designs |
