#!/usr/bin/env python3
"""L2 Test: Integration and cross-references between files."""
import os
import sys
import yaml
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Cross-File Integration")

    # 1. Decision matrix tool references exist on disk
    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    tool_refs = set()
    for section in ["peak_calling", "diff_binding", "annotation"]:
        for rule in dm.get(section, []):
            if "tool" in rule:
                tool_refs.add(rule["tool"])

    result.assert_true(len(tool_refs) >= 2, f"At least 2 tool references found (got {len(tool_refs)})")
    for ref in tool_refs:
        full_path = os.path.join(SKILL_DIR, ref)
        result.assert_file_exists(full_path, f"Tool ref '{ref}' exists")

    # 2. SKILL.md mentions all catalog files
    skill_path = os.path.join(SKILL_DIR, "SKILL.md")
    with open(skill_path, "r") as f:
        skill = f.read()

    catalog_files = [
        "preprocessing.md", "peak-calling.md", "diff-binding.md",
        "motif-enrichment.md", "annotation.md", "visualization.md"
    ]
    for cf in catalog_files:
        result.assert_true(f"tool-catalog/{cf}" in skill, f"SKILL.md references tool-catalog/{cf}")

    # 3. SKILL.md references decision-matrix.yaml
    result.assert_true("decision-matrix.yaml" in skill, "SKILL.md references decision-matrix.yaml")

    # 4. analyst-notebook.md references key tools
    nb_path = os.path.join(SKILL_DIR, "analyst-notebook.md")
    with open(nb_path, "r") as f:
        nb = f.read()
    for tool in ["MACS2", "DiffBind", "Bowtie2", "ChIPseeker", "FastQC", "MEME"]:
        result.assert_true(tool in nb, f"Notebook references {tool}")

    # 5. References file mentions key databases
    ref_path = os.path.join(SKILL_DIR, "references", "chipseq-plant-special.md")
    with open(ref_path, "r") as f:
        ref = f.read()
    for term in ["MACS2", "拟南芥", "抗体", "组蛋白"]:
        result.assert_true(term in ref, f"References mentions '{term}'")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
