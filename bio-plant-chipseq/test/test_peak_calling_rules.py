#!/usr/bin/env python3
"""L2 Test: Peak calling rule conditions - TF vs Histone."""
import os
import sys
import yaml
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Peak Calling Rule Conditions")

    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    peak_rules = dm["peak_calling"]

    # Both rules present
    r_tf = [r for r in peak_rules if r["rule_id"] == "chip-macs2-tf-001"][0]
    r_hist = [r for r in peak_rules if r["rule_id"] == "chip-macs2-histone-002"][0]

    # Conditions are distinct
    result.assert_equal(r_tf["condition"]["target_type"], "transcription_factor",
                        "TF rule targets transcription_factor")
    result.assert_equal(r_hist["condition"]["target_type"], "histone_mark",
                        "Histone rule targets histone_mark")
    result.assert_true(
        r_tf["condition"]["target_type"] != r_hist["condition"]["target_type"],
        "TF and histone conditions are distinct"
    )

    # Recommendations are distinct
    result.assert_equal(r_tf["recommend"], "macs2_narrow_peak",
                        "TF recommends MACS2 narrow peak")
    result.assert_equal(r_hist["recommend"], "macs2_broad_peak",
                        "Histone recommends MACS2 broad peak")
    result.assert_true(
        r_tf["recommend"] != r_hist["recommend"],
        "Narrow and broad peak recommendations are distinct"
    )

    # Both have reasons that are meaningful
    result.assert_true("narrow" in r_tf["reason"].lower() or "TF" in r_tf["reason"],
                       f"TF rule reason is meaningful: {r_tf['reason']}")
    result.assert_true("broad" in r_hist["reason"].lower() or "组蛋白" in r_hist["reason"],
                       f"Histone rule reason is meaningful: {r_hist['reason']}")

    # Both point to the same tool catalog entry
    result.assert_equal(r_tf["tool"], "tool-catalog/peak-calling.md", "TF rule points to peak-calling.md")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
