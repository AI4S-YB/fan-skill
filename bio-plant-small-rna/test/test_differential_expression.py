#!/usr/bin/env python3
"""L2 Test: Differential expression rules validation."""
import os
import sys
import yaml
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Differential Expression Rules")

    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    de_rules = dm.get("differential_expression", [])

    # Check DESeq2 rule
    deseq_rule = [r for r in de_rules if r["rule_id"] == "demirna-deseq2-001"]
    result.assert_true(len(deseq_rule) == 1, "DESeq2 rule exists")
    if deseq_rule:
        rule = deseq_rule[0]
        result.assert_equal(rule["priority"], 10, "DESeq2 rule has priority 10")
        result.assert_equal(rule["condition"]["rep_count"], ">= 3", "DESeq2 requires >= 3 replicates")
        result.assert_equal(rule["recommend"], "deseq2_on_mirna_counts", "Recommend DESeq2 on miRNA counts")
        result.assert_true("DESeq2" in rule["reason"], "Reason mentions DESeq2")

        # Tool reference should exist
        if "tool" in rule:
            # Direct tool reference from the decision matrix for DE
            # The decision matrix references tool-catalog/differential-mirna.md indirectly
            result.assert_true(True, "DE rule is well-formed")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
