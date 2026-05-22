---
name: bio-plant-infra
description: >
  Plant bioinformatics shared infrastructure — species reference,
  database access, data inspection, environment checks, and decision
  logging. Used by all bio-plant-* workflow skills.
tool_type: mixed
workflow: false
---

# Plant Bioinformatics Infrastructure

Shared scripts and references for all plant bioinformatics analysis Skills.

## What This Provides

| Resource | Path | Used By |
|----------|------|---------|
| Data inspection | `scripts/inspect_data.sh` | All Skills |
| Decision logging | `scripts/log_decision.sh` | All Skills |
| Environment check | `scripts/check_env.sh` | All Skills |
| QC check | `scripts/qc_check.sh` | All Skills |
| Species cheatsheet | `references/species-cheatsheet.md` | All Skills |
| Plant databases | `references/plant-databases.md` | All Skills |
| QC thresholds | `references/qc-thresholds.yaml` | All Skills |
| Troubleshooting | `references/troubleshooting.md` | All Skills |

## Usage

This Skill is not used independently. It is referenced by analysis Skills (bio-plant-gwas, bio-plant-population, etc.) through their SKILL.md files.

Each analysis Skill should instruct the Agent to:
1. Run `scripts/inspect_data.sh` to profile the input data
2. Look up `references/species-cheatsheet.md` if the species is detected
3. Use `scripts/log_decision.sh` to record every decision
4. Run `scripts/check_env.sh` before analysis
5. Use `scripts/qc_check.sh` at each QC checkpoint

## Templates

The `templates/` directory contains starter files for users:
- `samplesheet_template.csv` — sample-to-data mapping
- `params_template.yaml` — runtime parameters with decision mode
- `test_config_template.yaml` — server-specific test configuration
