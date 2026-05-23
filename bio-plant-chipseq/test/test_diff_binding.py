#!/usr/bin/env python3
"""L2 Test: Differential binding rules."""
import os
import sys
import yaml
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Differential Binding Rules")

    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    db_rules = dm.get("diff_binding", [])
    result.assert_true(len(db_rules) >= 1, f"Has diff binding rules (got {len(db_rules)})")

    db_rule = db_rules[0]
    result.assert_equal(db_rule["rule_id"], "chip-diffbind-001", "Rule ID is chip-diffbind-001")
    result.assert_equal(db_rule["priority"], 10, "Priority is 10")
    result.assert_equal(db_rule["condition"]["rep_count"], ">= 2", "Requires >= 2 replicates")
    result.assert_equal(db_rule["recommend"], "diffbind_deseq2", "Recommends DiffBind+DESeq2")
    result.assert_true("DiffBind" in db_rule["reason"], "Reason mentions DiffBind")
    result.assert_equal(db_rule["tool"], "tool-catalog/diff-binding.md", "Points to diff-binding.md")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
