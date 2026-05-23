# CRISPR Validation Strategy

**Goal:** Design and execute genotyping assays to confirm CRISPR/Cas9 editing events in plants
**Best for:** Characterizing edits in T0, T1, and T2 generations and screening for transgene-free mutants

## Prerequisites

- Target locus sequence (FASTA, including 500 bp up/downstream of cut site)
- Plant genomic DNA from edited plants
- Sanger sequencing or NGS capability

## Validation Workflow by Generation

### T0 Generation (Primary Transformants)

```python
# T0 plants are typically chimeric — multiple genotypes per plant
# Method: PCR + Sanger + ICE analysis

validation_pipeline_T0 = {
    "step1": "PCR amplify target locus (300-500 bp amplicon, cut site centered)",
    "step2": "Sanger sequence PCR product",
    "step3": "ICE (Inference of CRISPR Edits) analysis to estimate editing efficiency",
    "step4": "If ICE indicates edits: clone PCR product and sequence 10-20 colonies",
    "step5": "Report: editing types, frequencies, chimerism level"
}
```

### T1 Generation (Segregation)

```python
# Screen for homozygous/biallelic mutants and transgene-free plants
validation_pipeline_T1 = {
    "step1": "PCR amplify target locus for all T1 plants",
    "step2": "Sanger sequence (or T7E1 assay for quick screen)",
    "step3": "Classify: WT, heterozygous, homozygous, biallelic",
    "step4": "PCR for Cas9/sgRNA transgene (select transgene-free mutants)",
    "step5": "Select transgene-free homozygous mutants for T2 advancement"
}
```

## Genotyping Methods

### T7E1 / Surveyor Assay (Quick Screen)

```r
# T7E1 mismatch cleavage assay
# Semi-quantitative, good for screening large populations

# Expected bands:
# WT: single band (full-length PCR)
# Heterozygous: 3 bands (WT + 2 cleavage products)
# Homozygous/biallelic: 2 bands (2 cleavage products, no WT)
```

### ICE Analysis (Sanger Deconvolution)

```python
# Using Synthego ICE tool (https://ice.synthego.com/)
# or local analysis:

def run_ice_analysis(control_ab1, edited_ab1):
    """
    ICE: Inference of CRISPR Edits

    Input: Sanger .ab1 files for control (WT) and edited sample
    Output: Estimated indel distribution and editing efficiency
    """

    # Submit to ICE API or use local implementation
    # ice_results = ice_api.analyze(control_ab1, edited_ab1)

    return {
        "editing_efficiency": 0.75,  # % of reads with indels
        "indel_distribution": {
            "+1": 0.30,   # 1bp insertion
            "-1": 0.15,   # 1bp deletion
            "-3": 0.10,   # 3bp deletion
            "-7": 0.10,   # etc.
            "other": 0.10,
        },
        "ko_score": 0.65,  # Frameshift probability
    }
```

## PCR Primer Design for Genotyping

```python
def design_genotyping_primers(target_locus, cut_site_pos, product_size=400):
    """
    Design primers for amplifying CRISPR target site
    """
    # Amplicon should be centered on cut site
    flank = product_size // 2

    left_start = max(0, cut_site_pos - flank)
    right_end = min(len(target_locus), cut_site_pos + flank)

    # Design primers 50-100bp away from cut site
    primer_left_start = max(0, cut_site_pos - flank + 50)
    primer_right_start = min(len(target_locus), cut_site_pos + flank - 50)

    # Primer design (using primer3 or similar)
    from primer3 import design_primers

    left_seq = target_locus[primer_left_start:primer_left_start+200]
    right_seq = target_locus[primer_right_start-200:primer_right_start]

    primers = design_primers(
        seq_template=target_locus,
        target_start=primer_left_start,
        target_end=primer_right_start,
        product_size_range=(product_size - 50, product_size + 50)
    )

    return primers
```

## Transgene Detection

```python
# PCR for Cas9 and selectable marker to identify transgene-free plants
transgene_primers = {
    "SpCas9": {
        "forward": "GACAAGAAGTACAGCATCGGCC",
        "reverse": "GTTGATGATGTAGTCGCTGGTG",
        "product_size": 450
    },
    "Hygromycin": {
        "forward": "GATGTTGGCGACCTCGTATT",
        "reverse": "CCATGTAGTGTATTGACCGATTCC",
        "product_size": 600
    },
    "BAR": {
        "forward": "TCTGCACCATCGTCAACCACTA",
        "reverse": "TAGATCTTCGGTGACGGGCAG",
        "product_size": 350
    }
}

# Also include endogenous control PCR (e.g., Actin, Ubiquitin)
# If transgene PCR = negative, and control PCR = positive
# -> transgene-free plant
```

## Key Parameters

| Parameter | Recommended | Rationale |
|-----------|------------|-----------|
| PCR amplicon size | 300-500 bp | Clear separation of WT vs edited bands |
| Distance of primers from cut | > 50 bp | Avoid large deletions removing primer sites |
| ICE R^2 | > 0.8 | Good fit indicates reliable estimation |
| Cloning colonies (T0) | 10-20 | Capture chimeric diversity |
| T1 screening | 20-50 plants | Mendelian segregation |
| Transgene-free screening | T1 or T2 | Segregation needed |

## Plant-Specific Notes

- T0 genotyping can be misleading: a plant may be chimeric, showing edits in one leaf but not another
- Sample T0 plants from 2-3 different leaves for chimerism assessment
- T1 segregation: screen >= 20 plants to find transgene-free homozygous mutants
- For vegetatively propagated crops (potato, cassava): once an edited line is established, it stays (no segregation)
- For seed-propagated crops: T2 = first generation of true-breeding homozygous edited line
- Polyploid editing validation: PCR + cloning + sequencing 50+ colonies to capture all subgenome combinations

## Comprehensive Validation Report

```python
def generate_validation_report(plant_id, generation, genotyping_data):
    """Generate CRISPR editing validation report"""
    report = []

    report.append(f"=== CRISPR Validation Report ===")
    report.append(f"Plant ID: {plant_id}")
    report.append(f"Generation: {generation}")
    report.append(f"")

    report.append(f"Target gene: {genotyping_data['target_gene']}")
    report.append(f"sgRNA: {genotyping_data['sgrna_sequence']}")
    report.append(f"")

    report.append(f"Editing efficiency: {genotyping_data['editing_efficiency']:.1%}")
    report.append(f"Frameshift probability: {genotyping_data['ko_score']:.1%}")
    report.append(f"")

    report.append("Indel Distribution:")
    for indel, freq in genotyping_data['indel_distribution'].items():
        report.append(f"  {indel}: {freq:.1%}")

    report.append(f"")
    report.append(f"Genotype: {genotyping_data['genotype']}")
    report.append(f"Transgene status: {genotyping_data['transgene_status']}")

    return "\n".join(report)
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| PCR fails | Large deletion removed primer site | Use primers further away from cut site |
| ICE R^2 < 0.5 | Poor sequencing quality or complex editing | Clone PCR and sequence colonies |
| T7E1 shows cleavage in WT | Natural indel polymorphism at locus | Sequence WT first to confirm clean target site |
| "Wild-type" T1 plants >75% | Chimeric T0 parent, low germline editing | Screen more T1 plants, use T0 plants with higher editing |
| Cannot find transgene-free mutants | Cas9 and sgRNA on same T-DNA | Use separate T-DNA design for marker and editing components |
