#!/usr/bin/env python3
"""
L2 Test Engine for bio-plant-chipseq Skill.
Tests structure, metadata, decision matrix rules, tool catalog completeness,
references, analyst notebook, and integration.
"""

import os
import sys
import yaml
import re

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
    required_dirs = ["tool-catalog", "references", "test"]
    for d in required_dirs:
        result.assert_file_exists(os.path.join(SKILL_DIR, d), f"Directory '{d}' exists")

    required_files = ["SKILL.md", "decision-matrix.yaml", "analyst-notebook.md"]
    for f in required_files:
        result.assert_file_exists(os.path.join(SKILL_DIR, f), f"File '{f}' exists")


def test_skill_metadata(result):
    """Test SKILL.md metadata."""
    print("\n[2] SKILL.md Metadata Tests")
    skill_path = os.path.join(SKILL_DIR, "SKILL.md")
    with open(skill_path, "r") as f:
        content = f.read()

    result.assert_true("skill_name: bio-plant-chipseq" in content, "Has skill_name field")
    result.assert_true("version:" in content, "Has version field")
    result.assert_true("category: plant-bioinformatics" in content, "Has category field")
    result.assert_true("analysis_type: chip_sequencing" in content, "Has analysis_type field")
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
    result.assert_equal(dm.get("analysis"), "chipseq", "Analysis type is chipseq")

    # Test peak calling rules
    peak_rules = dm.get("peak_calling", [])
    result.assert_true(len(peak_rules) >= 2, f"Has at least 2 peak calling rules (got {len(peak_rules)})")

    rule_ids = [r["rule_id"] for r in peak_rules]
    result.assert_true("chip-macs2-tf-001" in rule_ids, "Has TF narrow peak rule")
    result.assert_true("chip-macs2-histone-002" in rule_ids, "Has histone broad peak rule")

    # Check TF rule
    tf_rule = [r for r in peak_rules if r["rule_id"] == "chip-macs2-tf-001"][0]
    result.assert_equal(tf_rule["priority"], 10, "TF rule priority is 10")
    result.assert_equal(tf_rule["condition"]["target_type"], "transcription_factor", "TF condition is transcription_factor")
    result.assert_equal(tf_rule["recommend"], "macs2_narrow_peak", "TF recommends narrow peak")
    result.assert_true("tool" in tf_rule, "TF rule has tool reference")

    # Check histone rule
    hist_rule = [r for r in peak_rules if r["rule_id"] == "chip-macs2-histone-002"][0]
    result.assert_equal(hist_rule["priority"], 10, "Histone rule priority is 10")
    result.assert_equal(hist_rule["condition"]["target_type"], "histone_mark", "Histone condition is histone_mark")
    result.assert_equal(hist_rule["recommend"], "macs2_broad_peak", "Histone recommends broad peak")

    # Test diff binding rules
    db_rules = dm.get("diff_binding", [])
    result.assert_true(len(db_rules) >= 1, f"Has at least 1 diff binding rule (got {len(db_rules)})")
    db_rule = db_rules[0]
    result.assert_equal(db_rule["rule_id"], "chip-diffbind-001", "Diff bind rule ID correct")
    result.assert_equal(db_rule["condition"]["rep_count"], ">= 2", "Diff bind requires >= 2 reps")

    # Test annotation rules
    anno_rules = dm.get("annotation", [])
    result.assert_true(len(anno_rules) >= 1, f"Has at least 1 annotation rule (got {len(anno_rules)})")
    anno_rule = anno_rules[0]
    result.assert_equal(anno_rule["rule_id"], "chip-annotate-plant-001", "Annotation rule ID correct")
    result.assert_equal(anno_rule["priority"], 5, "Annotation rule priority is 5")

    # Test fallback
    fallback = dm.get("fallback", [])
    result.assert_true(len(fallback) >= 1, f"Has at least 1 fallback rule (got {len(fallback)})")
    fb_rule = fallback[0]
    result.assert_equal(fb_rule["rule_id"], "chip-fallback-999", "Fallback rule ID correct")
    result.assert_equal(fb_rule["priority"], 0, "Fallback priority is 0")
    result.assert_equal(fb_rule.get("action"), "delegate_to_expert", "Fallback action is delegate_to_expert")


def test_decision_matrix_conditions(result):
    """Test specific conditions in decision matrix."""
    print("\n[4] Decision Matrix Conditions Tests")
    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    peak_rules = dm.get("peak_calling", [])
    tf_rule = [r for r in peak_rules if r["rule_id"] == "chip-macs2-tf-001"][0]
    hist_rule = [r for r in peak_rules if r["rule_id"] == "chip-macs2-histone-002"][0]

    result.assert_equal(tf_rule["condition"]["target_type"], "transcription_factor",
                        "TF rule has transcription_factor condition")
    result.assert_equal(hist_rule["condition"]["target_type"], "histone_mark",
                        "Histone rule has histone_mark condition")

    # Both peak calling rules have equal priority (both valid, user chooses)
    result.assert_equal(tf_rule["priority"], hist_rule["priority"],
                        "TF and histone rules have same priority")

    # Diff bind condition
    db_rule = dm["diff_binding"][0]
    result.assert_equal(db_rule["condition"]["rep_count"], ">= 2",
                        "Diff bind requires >= 2 replicates")


def test_tool_catalog(result):
    """Test tool catalog files exist and have content."""
    print("\n[5] Tool Catalog Tests")
    catalog_dir = os.path.join(SKILL_DIR, "tool-catalog")

    expected_files = [
        "peak-calling.md",
        "diff-binding.md",
        "motif-enrichment.md",
        "annotation.md",
        "preprocessing.md",
        "visualization.md",
    ]

    for f in expected_files:
        path = os.path.join(catalog_dir, f)
        result.assert_file_exists(path, f"Tool catalog file '{f}' exists")
        result.assert_file_not_empty(path, f"Tool catalog file '{f}' is not empty")

        with open(path, "r") as fh:
            content = fh.read()
        line_count = len(content.strip().split("\n"))
        result.assert_true(line_count >= 10, f"'{f}' has at least 10 lines of content (got {line_count})")
        result.assert_true("#" in content, f"'{f}' has markdown headings")


def test_tool_catalog_content_quality(result):
    """Test that tool catalog files contain specific expected content."""
    print("\n[6] Tool Catalog Content Quality Tests")
    catalog_dir = os.path.join(SKILL_DIR, "tool-catalog")

    # Check preprocessing.md
    with open(os.path.join(catalog_dir, "preprocessing.md"), "r") as f:
        preproc = f.read()
    result.assert_true("FastQC" in preproc, "preprocessing.md mentions FastQC")
    result.assert_true("Bowtie2" in preproc, "preprocessing.md mentions Bowtie2")
    result.assert_true("SAMtools" in preproc, "preprocessing.md mentions SAMtools")
    result.assert_true("植物" in preproc, "preprocessing.md has plant-specific content")

    # Check peak-calling.md
    with open(os.path.join(catalog_dir, "peak-calling.md"), "r") as f:
        peak = f.read()
    result.assert_true("MACS2" in peak, "peak-calling.md mentions MACS2")
    result.assert_true("narrow" in peak.lower(), "peak-calling.md mentions narrow peak")
    result.assert_true("broad" in peak.lower(), "peak-calling.md mentions broad peak")
    result.assert_true("transcription_factor" in peak.lower() or "转录因子" in peak,
                       "peak-calling.md discusses TF peaks")

    # Check diff-binding.md
    with open(os.path.join(catalog_dir, "diff-binding.md"), "r") as f:
        db = f.read()
    result.assert_true("DiffBind" in db, "diff-binding.md mentions DiffBind")
    result.assert_true("DESeq2" in db, "diff-binding.md mentions DESeq2")

    # Check motif-enrichment.md
    with open(os.path.join(catalog_dir, "motif-enrichment.md"), "r") as f:
        me = f.read()
    result.assert_true("MEME" in me, "motif-enrichment.md mentions MEME")
    result.assert_true("motif" in me.lower(), "motif-enrichment.md discusses motifs")

    # Check annotation.md
    with open(os.path.join(catalog_dir, "annotation.md"), "r") as f:
        anno = f.read()
    result.assert_true("ChIPseeker" in anno, "annotation.md mentions ChIPseeker")

    # Check visualization.md
    with open(os.path.join(catalog_dir, "visualization.md"), "r") as f:
        viz = f.read()
    result.assert_true("IGV" in viz or "deepTools" in viz, "visualization.md mentions visualization tools")


def test_references(result):
    """Test references directory."""
    print("\n[7] References Tests")
    ref_path = os.path.join(SKILL_DIR, "references", "chipseq-plant-special.md")
    result.assert_file_exists(ref_path, "references/chipseq-plant-special.md exists")
    result.assert_file_not_empty(ref_path, "references/chipseq-plant-special.md is not empty")

    with open(ref_path, "r") as f:
        ref_content = f.read()

    result.assert_true("植物" in ref_content, "References contains Chinese content")
    result.assert_true("MACS2" in ref_content or "macs2" in ref_content.lower(), "References mentions MACS2")
    result.assert_true("组蛋白" in ref_content or "histone" in ref_content.lower(), "References mentions histone")
    result.assert_true("抗体" in ref_content or "antibody" in ref_content.lower(), "References mentions antibody")
    result.assert_true("转座子" in ref_content or "TE" in ref_content, "References mentions transposons")
    result.assert_true("拟南芥" in ref_content, "References mentions Arabidopsis")

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

    result.assert_true("植物" in nb_content, "Notebook contains Chinese content")
    result.assert_true("ChIP-seq" in nb_content or "ChIP" in nb_content, "Notebook covers ChIP-seq")
    result.assert_true("MACS2" in nb_content, "Notebook mentions MACS2")
    result.assert_true("DiffBind" in nb_content, "Notebook mentions DiffBind")
    result.assert_true("Bowtie2" in nb_content, "Notebook mentions Bowtie2")
    result.assert_true("```" in nb_content, "Notebook contains code blocks")

    sections = ["预处理", "Peak Calling", "Peak 注释", "Motif", "差异", "可视化"]
    found = 0
    for s in sections:
        if s.lower() in nb_content.lower():
            found += 1
    result.assert_true(found >= 4, f"Notebook covers major analysis sections (found {found}/6)")

    line_count = len(nb_content.strip().split("\n"))
    result.assert_true(line_count >= 80, f"Notebook has at least 80 lines (got {line_count})")


def test_integration(result):
    """Integration test: cross-reference between files."""
    print("\n[9] Integration Tests")

    # Check that decision matrix tool references point to existing files
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
        ref_path = os.path.join(SKILL_DIR, ref)
        result.assert_file_exists(ref_path, f"Decision matrix tool ref '{ref}' points to existing file")

    # Check SKILL.md tool catalog table references match actual files
    skill_path = os.path.join(SKILL_DIR, "SKILL.md")
    with open(skill_path, "r") as f:
        skill_content = f.read()

    catalog_files = [
        "preprocessing.md",
        "peak-calling.md",
        "diff-binding.md",
        "motif-enrichment.md",
        "annotation.md",
        "visualization.md",
    ]
    for cf in catalog_files:
        result.assert_true(cf in skill_content, f"SKILL.md references tool-catalog/{cf}")

    # Check notebook references key tools
    nb_path = os.path.join(SKILL_DIR, "analyst-notebook.md")
    with open(nb_path, "r") as f:
        nb = f.read()
    for tool in ["MACS2", "DiffBind", "Bowtie2", "ChIPseeker"]:
        result.assert_true(tool in nb, f"Notebook references {tool}")


def test_rule_priority_ordering(result):
    """Test that rules within each category are ordered by priority."""
    print("\n[10] Rule Priority Ordering Tests")
    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    for section in ["peak_calling", "diff_binding", "annotation"]:
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
        for section in ["peak_calling", "diff_binding", "annotation"]:
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
