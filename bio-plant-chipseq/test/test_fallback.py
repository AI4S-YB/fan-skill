#!/usr/bin/env python3
"""L2 Test: Fallback rule validation."""
import os
import sys
import yaml
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Fallback Rule Validation")

    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    fallback = dm.get("fallback", [])
    result.assert_equal(len(fallback), 1, "Exactly 1 fallback rule")

    fb = fallback[0]
    result.assert_equal(fb["rule_id"], "chip-fallback-999", "Fallback rule ID is chip-fallback-999")
    result.assert_equal(fb["priority"], 0, "Fallback priority is 0")
    result.assert_equal(fb.get("condition"), "none_matched", "Fallback condition is 'none_matched'")
    result.assert_equal(fb.get("action"), "delegate_to_expert", "Fallback action is 'delegate_to_expert'")

    # Fallback should be lowest priority across all rules
    all_priorities = []
    for section in ["peak_calling", "diff_binding", "annotation"]:
        for rule in dm.get(section, []):
            all_priorities.append(rule["priority"])
    min_rule_priority = min(all_priorities) if all_priorities else 0
    result.assert_true(fb["priority"] < min_rule_priority,
                       f"Fallback priority ({fb['priority']}) < minimum rule priority ({min_rule_priority})")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
