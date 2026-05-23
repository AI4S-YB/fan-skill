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

    result.assert_equal(dm.get("version"), "1.0", "Version is 1.0")
    result.assert_equal(dm.get("analysis"), "chipseq", "Analysis type is chipseq")

    # Required sections
    for section in ["peak_calling", "diff_binding", "annotation", "fallback"]:
        result.assert_true(section in dm, f"Has '{section}' section")
        result.assert_true(len(dm[section]) >= 1, f"'{section}' has at least 1 rule")

    # Peak calling rules
    peak_rules = dm["peak_calling"]
    rule_ids = [r["rule_id"] for r in peak_rules]
    result.assert_true("chip-macs2-tf-001" in rule_ids, "Has TF rule")
    result.assert_true("chip-macs2-histone-002" in rule_ids, "Has histone rule")

    # TF rule
    r_tf = [r for r in peak_rules if r["rule_id"] == "chip-macs2-tf-001"][0]
    for field in ["rule_id", "priority", "condition", "recommend", "reason", "tool"]:
        result.assert_true(field in r_tf, f"TF rule has '{field}' field")
    result.assert_equal(r_tf["condition"]["target_type"], "transcription_factor", "TF condition correct")
    result.assert_equal(r_tf["recommend"], "macs2_narrow_peak", "TF recommends narrow peak")

    # Histone rule
    r_hist = [r for r in peak_rules if r["rule_id"] == "chip-macs2-histone-002"][0]
    result.assert_equal(r_hist["condition"]["target_type"], "histone_mark", "Histone condition correct")
    result.assert_equal(r_hist["recommend"], "macs2_broad_peak", "Histone recommends broad peak")

    # Diff bind rule
    r_db = dm["diff_binding"][0]
    result.assert_equal(r_db["rule_id"], "chip-diffbind-001", "Diff bind rule ID correct")
    result.assert_equal(r_db["priority"], 10, "Diff bind priority is 10")

    # Annotation rule
    r_anno = dm["annotation"][0]
    result.assert_equal(r_anno["rule_id"], "chip-annotate-plant-001", "Annotation rule ID correct")
    result.assert_equal(r_anno["priority"], 5, "Annotation priority is 5")

    # Fallback
    fb = dm["fallback"][0]
    result.assert_equal(fb["rule_id"], "chip-fallback-999", "Fallback rule ID correct")
    result.assert_equal(fb["priority"], 0, "Fallback priority is 0")
    result.assert_equal(fb["action"], "delegate_to_expert", "Fallback action is delegate_to_expert")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
