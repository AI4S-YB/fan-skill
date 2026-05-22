# Plant Marker Development Specifics

## Marker Types by Breeding Application

| Application | Throughput | Best Marker |
|------------|:---------:|------------|
| Foreground selection (1-5 genes) | Low | KASP or gel-based |
| Background selection (genome-wide) | Medium | SNP chip or GBS |
| Quality assurance (seed purity) | Low | SSR or InDel |
| Diagnostic marker (specific allele) | Low | KASP or dCAPS |

## Polyploid Marker Design

### Allopolyploids
- Design subgenome-specific primers using subgenome-diagnostic SNPs
- Validate by BLAST against each subgenome reference

### Autopolyploids
- No single-copy design possible — use dosage scoring
- KASP can still work with allele-specific fluorescence ratio
