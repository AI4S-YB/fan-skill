#!/usr/bin/env python3
"""L2 Test: Decision matrix rules and conditions."""
import os
import sys
import yaml
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Decision Matrix Rules")

    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    # Version and type
    result.assert_equal(dm.get("version"), "1.0", "Version is 1.0")
    result.assert_equal(dm.get("analysis"), "small_rna", "Analysis type is small_rna")

    # Required sections exist
    for section in ["mirna_prediction", "target_prediction", "differential_expression", "fallback"]:
        result.assert_true(section in dm, f"Has '{section}' section")
        result.assert_true(len(dm[section]) >= 1, f"'{section}' has at least 1 rule")

    # miRNA prediction rules
    mirna_rules = dm["mirna_prediction"]
    rule_ids = [r["rule_id"] for r in mirna_rules]
    result.assert_true("mirna-plant-species-001" in rule_ids, "Has model species rule")
    result.assert_true("mirna-novel-002" in rule_ids, "Has novel discovery rule")

    # Verify rule 001 structure
    r001 = [r for r in mirna_rules if r["rule_id"] == "mirna-plant-species-001"][0]
    for field in ["rule_id", "priority", "condition", "recommend", "reason", "tool"]:
        result.assert_true(field in r001, f"Rule 001 has '{field}' field")

    # Verify rule priorities
    result.assert_true(r001["priority"] > 0, "Model species rule has positive priority")

    # Target prediction
    target_rule = dm["target_prediction"][0]
    result.assert_equal(target_rule["rule_id"], "target-psrna-001", "Target rule ID correct")
    result.assert_equal(target_rule["priority"], 10, "Target rule priority is 10")

    # Differential expression
    de_rule = dm["differential_expression"][0]
    result.assert_equal(de_rule["rule_id"], "demirna-deseq2-001", "DE rule ID correct")

    # Fallback
    fb = dm["fallback"][0]
    result.assert_equal(fb["rule_id"], "smrna-fallback-999", "Fallback rule ID correct")
    result.assert_equal(fb["priority"], 0, "Fallback priority is 0")
    result.assert_equal(fb["action"], "delegate_to_expert", "Fallback action is delegate_to_expert")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
