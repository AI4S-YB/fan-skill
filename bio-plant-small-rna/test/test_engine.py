#!/usr/bin/env python3
"""
L2 Test Engine for bio-plant-small-rna Skill.
Tests structure, metadata, decision matrix rules, tool catalog completeness,
references, analyst notebook, and integration.
"""

import os
import sys
import yaml
import re
import json

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class TestResult:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.errors = []

    def assert_true(self, condition, message):
        if condition:
            self.passed += 1
            print(f"  PASS: {message}")
        else:
            self.failed += 1
            self.errors.append(message)
            print(f"  FAIL: {message}")

    def assert_equal(self, actual, expected, message):
        if actual == expected:
            self.passed += 1
            print(f"  PASS: {message}")
        else:
            self.failed += 1
            err = f"{message} (expected={expected}, got={actual})"
            self.errors.append(err)
            print(f"  FAIL: {err}")

    def assert_file_exists(self, path, message):
        if os.path.exists(path):
            self.passed += 1
            print(f"  PASS: {message}")
        else:
            self.failed += 1
            self.errors.append(message)
            print(f"  FAIL: {message}")

    def assert_file_not_empty(self, path, message):
        if os.path.exists(path) and os.path.getsize(path) > 0:
            self.passed += 1
            print(f"  PASS: {message}")
        else:
            self.failed += 1
            self.errors.append(message)
            print(f"  FAIL: {message}")

    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'='*60}")
        print(f"TEST SUMMARY: {self.passed}/{total} passed, {self.failed} failed")
        if self.errors:
            print(f"\nFailed tests:")
            for e in self.errors:
                print(f"  - {e}")
        print(f"{'='*60}")
        return self.failed == 0


def test_structure(result):
    """Test directory and file structure."""
    print("\n[1] Structure Tests")
    required_dirs = [
        "tool-catalog",
        "references",
        "test",
    ]
    for d in required_dirs:
        result.assert_file_exists(os.path.join(SKILL_DIR, d), f"Directory '{d}' exists")

    required_files = [
        "SKILL.md",
        "decision-matrix.yaml",
        "analyst-notebook.md",
    ]
    for f in required_files:
        result.assert_file_exists(os.path.join(SKILL_DIR, f), f"File '{f}' exists")


def test_skill_metadata(result):
    """Test SKILL.md metadata."""
    print("\n[2] SKILL.md Metadata Tests")
    skill_path = os.path.join(SKILL_DIR, "SKILL.md")
    with open(skill_path, "r") as f:
        content = f.read()

    result.assert_true("skill_name: bio-plant-small-rna" in content, "Has skill_name field")
    result.assert_true("version:" in content, "Has version field")
    result.assert_true("category: plant-bioinformatics" in content, "Has category field")
    result.assert_true("analysis_type: small_rna_sequencing" in content, "Has analysis_type field")
    result.assert_true("decision-matrix.yaml" in content, "References decision-matrix.yaml")
    result.assert_true("tool-catalog/" in content, "References tool-catalog/")
    result.assert_true("## Overview" in content or "## When to Use" in content, "Has structured sections")
    result.assert_true("## Decision Matrix" in content, "Has Decision Matrix section")
    result.assert_true("## Tool Catalog" in content, "Has Tool Catalog section")
    result.assert_true("## Fallback" in content, "Has Fallback section")


def test_decision_matrix_structure(result):
    """Test decision-matrix.yaml structure and rules."""
    print("\n[3] Decision Matrix Structure Tests")
    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    result.assert_equal(dm.get("version"), "1.0", "Version is 1.0")
    result.assert_equal(dm.get("analysis"), "small_rna", "Analysis type is small_rna")

    # Test miRNA prediction rules
    mirna_rules = dm.get("mirna_prediction", [])
    result.assert_true(len(mirna_rules) >= 2, f"Has at least 2 miRNA prediction rules (got {len(mirna_rules)})")

    rule_ids = [r["rule_id"] for r in mirna_rules]
    result.assert_true("mirna-plant-species-001" in rule_ids, "Has mirna-plant-species-001 rule")
    result.assert_true("mirna-novel-002" in rule_ids, "Has mirna-novel-002 rule")

    # Check rule 001 (model species)
    rule001 = [r for r in mirna_rules if r["rule_id"] == "mirna-plant-species-001"][0]
    result.assert_equal(rule001["priority"], 10, "mirna-plant-species-001 priority is 10")
    result.assert_true("recommend" in rule001, "Rule 001 has recommend field")
    result.assert_true("reason" in rule001, "Rule 001 has reason field")
    result.assert_true("tool" in rule001, "Rule 001 has tool field")

    # Check rule 002 (novel)
    rule002 = [r for r in mirna_rules if r["rule_id"] == "mirna-novel-002"][0]
    result.assert_equal(rule002["priority"], 8, "mirna-novel-002 priority is 8")

    # Test target prediction rules
    target_rules = dm.get("target_prediction", [])
    result.assert_true(len(target_rules) >= 1, f"Has at least 1 target prediction rule (got {len(target_rules)})")
    result.assert_true(any(r["rule_id"] == "target-psrna-001" for r in target_rules), "Has psRNATarget rule")

    # Test differential expression rules
    de_rules = dm.get("differential_expression", [])
    result.assert_true(len(de_rules) >= 1, f"Has at least 1 DE rule (got {len(de_rules)})")
    result.assert_true(any(r["rule_id"] == "demirna-deseq2-001" for r in de_rules), "Has DESeq2 DE rule")

    # Test fallback
    fallback = dm.get("fallback", [])
    result.assert_true(len(fallback) >= 1, f"Has at least 1 fallback rule (got {len(fallback)})")
    result.assert_true(any(r["rule_id"] == "smrna-fallback-999" for r in fallback), "Has fallback rule")
    fb_rule = [r for r in fallback if r["rule_id"] == "smrna-fallback-999"][0]
    result.assert_equal(fb_rule["priority"], 0, "Fallback priority is 0")
    result.assert_equal(fb_rule.get("action"), "delegate_to_expert", "Fallback action is delegate_to_expert")


def test_decision_matrix_conditions(result):
    """Test specific conditions in decision matrix."""
    print("\n[4] Decision Matrix Conditions Tests")
    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    mirna_rules = dm.get("mirna_prediction", [])
    rule001 = [r for r in mirna_rules if r["rule_id"] == "mirna-plant-species-001"][0]

    species = rule001["condition"].get("species_detected", [])
    result.assert_true("oryza_sativa" in species, "Species list includes oryza_sativa")
    result.assert_true("zea_mays" in species, "Species list includes zea_mays")
    result.assert_true("arabidopsis_thaliana" in species, "Species list includes arabidopsis_thaliana")
    result.assert_true("glycine_max" in species, "Species list includes glycine_max")
    result.assert_equal(len(species), 4, "Species list has exactly 4 model species")

    rule002 = [r for r in mirna_rules if r["rule_id"] == "mirna-novel-002"][0]
    result.assert_equal(rule002["condition"].get("species_detected"), "unknown", "Novel rule condition is 'unknown'")

    de_rules = dm.get("differential_expression", [])
    de_rule = [r for r in de_rules if r["rule_id"] == "demirna-deseq2-001"][0]
    result.assert_equal(de_rule["condition"].get("rep_count"), ">= 3", "DE rule requires >= 3 replicates")


def test_tool_catalog(result):
    """Test tool catalog files exist and have content."""
    print("\n[5] Tool Catalog Tests")
    catalog_dir = os.path.join(SKILL_DIR, "tool-catalog")

    expected_files = [
        "mirna-prediction.md",
        "target-prediction.md",
        "degradome.md",
        "differential-mirna.md",
        "preprocessing.md",
        "visualization.md",
    ]

    for f in expected_files:
        path = os.path.join(catalog_dir, f)
        result.assert_file_exists(path, f"Tool catalog file '{f}' exists")
        result.assert_file_not_empty(path, f"Tool catalog file '{f}' is not empty")

        # Check minimum content: should have at least a title and description
        with open(path, "r") as fh:
            content = fh.read()
        line_count = len(content.strip().split("\n"))
        result.assert_true(line_count >= 10, f"'{f}' has at least 10 lines of content (got {line_count})")
        result.assert_true("#" in content, f"'{f}' has markdown headings")


def test_tool_catalog_content_quality(result):
    """Test that tool catalog files contain specific expected content."""
    print("\n[6] Tool Catalog Content Quality Tests")

    catalog_dir = os.path.join(SKILL_DIR, "tool-catalog")

    # Check preprocessing.md mentions key tools
    with open(os.path.join(catalog_dir, "preprocessing.md"), "r") as f:
        preproc = f.read()
    result.assert_true("FastQC" in preproc, "preprocessing.md mentions FastQC")
    result.assert_true("Cutadapt" in preproc or "cutadapt" in preproc, "preprocessing.md mentions Cutadapt")
    result.assert_true("18" in preproc and "30" in preproc, "preprocessing.md mentions length range 18-30")

    # Check mirna-prediction.md
    with open(os.path.join(catalog_dir, "mirna-prediction.md"), "r") as f:
        mirna = f.read()
    result.assert_true("miRDeep2" in mirna or "mirdeep2" in mirna.lower(), "mirna-prediction.md mentions miRDeep2")

    # Check target-prediction.md
    with open(os.path.join(catalog_dir, "target-prediction.md"), "r") as f:
        target = f.read()
    result.assert_true("psRNATarget" in target, "target-prediction.md mentions psRNATarget")

    # Check degradome.md
    with open(os.path.join(catalog_dir, "degradome.md"), "r") as f:
        degradome = f.read()
    result.assert_true("CleaveLand" in degradome or "cleaveland" in degradome.lower(), "degradome.md mentions CleaveLand")

    # Check differential-mirna.md
    with open(os.path.join(catalog_dir, "differential-mirna.md"), "r") as f:
        demirna = f.read()
    result.assert_true("DESeq2" in demirna, "differential-mirna.md mentions DESeq2")

    # Check visualization.md
    with open(os.path.join(catalog_dir, "visualization.md"), "r") as f:
        viz = f.read()
    result.assert_true("ggplot2" in viz or "Volcano" in viz or "Heatmap" in viz, "visualization.md mentions key visualizations")


def test_references(result):
    """Test references directory."""
    print("\n[7] References Tests")
    ref_path = os.path.join(SKILL_DIR, "references", "small-rna-plant-special.md")
    result.assert_file_exists(ref_path, "references/small-rna-plant-special.md exists")
    result.assert_file_not_empty(ref_path, "references/small-rna-plant-special.md is not empty")

    with open(ref_path, "r") as f:
        ref_content = f.read()

    result.assert_true("miRNA" in ref_content or "microRNA" in ref_content, "References mentions miRNA")
    result.assert_true("植物" in ref_content, "References contains Chinese content")
    result.assert_true("siRNA" in ref_content, "References mentions siRNA")
    result.assert_true("miRBase" in ref_content, "References mentions miRBase")
    result.assert_true("降解组" in ref_content or "Degradome" in ref_content, "References mentions degradome")

    # Should have at least 50 lines
    line_count = len(ref_content.strip().split("\n"))
    result.assert_true(line_count >= 50, f"References has at least 50 lines (got {line_count})")


def test_analyst_notebook(result):
    """Test analyst notebook."""
    print("\n[8] Analyst Notebook Tests")
    nb_path = os.path.join(SKILL_DIR, "analyst-notebook.md")
    result.assert_file_exists(nb_path, "analyst-notebook.md exists")
    result.assert_file_not_empty(nb_path, "analyst-notebook.md is not empty")

    with open(nb_path, "r") as f:
        nb_content = f.read()

    result.assert_true("预处理" in nb_content or "preprocessing" in nb_content.lower(), "Notebook covers preprocessing")
    result.assert_true("miRNA" in nb_content or "microRNA" in nb_content, "Notebook covers miRNA analysis")
    result.assert_true("差异" in nb_content or "differential" in nb_content.lower(), "Notebook covers differential analysis")
    result.assert_true("DESeq2" in nb_content, "Notebook mentions DESeq2")
    result.assert_true("FastQC" in nb_content or "fastqc" in nb_content, "Notebook mentions FastQC")
    result.assert_true("psRNATarget" in nb_content, "Notebook mentions psRNATarget")

    line_count = len(nb_content.strip().split("\n"))
    result.assert_true(line_count >= 80, f"Notebook has at least 80 lines (got {line_count})")


def test_integration(result):
    """Integration test: cross-reference between files."""
    print("\n[9] Integration Tests")

    # Check that decision matrix tool references point to existing files
    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    # Collect all tool references from decision matrix
    tool_refs = []
    for section in ["mirna_prediction", "target_prediction", "differential_expression"]:
        for rule in dm.get(section, []):
            if "tool" in rule:
                tool_refs.append(rule["tool"])

    for ref in tool_refs:
        ref_path = os.path.join(SKILL_DIR, ref)
        result.assert_file_exists(ref_path, f"Decision matrix tool ref '{ref}' points to existing file")

    # Check SKILL.md tool catalog table references match actual files
    skill_path = os.path.join(SKILL_DIR, "SKILL.md")
    with open(skill_path, "r") as f:
        skill_content = f.read()

    catalog_files = [
        "preprocessing.md",
        "mirna-prediction.md",
        "target-prediction.md",
        "degradome.md",
        "differential-mirna.md",
        "visualization.md",
    ]
    for cf in catalog_files:
        result.assert_true(cf in skill_content, f"SKILL.md references tool-catalog/{cf}")


def test_rule_priority_ordering(result):
    """Test that rules within each category are ordered by priority."""
    print("\n[10] Rule Priority Ordering Tests")
    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    for section in ["mirna_prediction", "target_prediction", "differential_expression"]:
        rules = dm.get(section, [])
        if len(rules) > 1:
            for i in range(len(rules) - 1):
                result.assert_true(
                    rules[i]["priority"] >= rules[i+1]["priority"],
                    f"{section} rule {rules[i]['rule_id']} priority ({rules[i]['priority']}) >= "
                    f"{rules[i+1]['rule_id']} priority ({rules[i+1]['priority']})"
                )

    # Fallback should have lowest priority
    fallback = dm.get("fallback", [])
    if fallback:
        fb_pri = fallback[0]["priority"]
        for section in ["mirna_prediction", "target_prediction", "differential_expression"]:
            for rule in dm.get(section, []):
                result.assert_true(
                    rule["priority"] > fb_pri,
                    f"Rule {rule['rule_id']} priority ({rule['priority']}) > fallback priority ({fb_pri})"
                )


def main():
    print(f"Testing Skill: {SKILL_DIR}")
    print(f"{'='*60}")

    result = TestResult()

    test_structure(result)
    test_skill_metadata(result)
    test_decision_matrix_structure(result)
    test_decision_matrix_conditions(result)
    test_tool_catalog(result)
    test_tool_catalog_content_quality(result)
    test_references(result)
    test_analyst_notebook(result)
    test_integration(result)
    test_rule_priority_ordering(result)

    return result.summary()


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
