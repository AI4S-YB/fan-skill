#!/usr/bin/env python3
"""L2 Test: References files validation."""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: References Validation")

    ref_path = os.path.join(SKILL_DIR, "references", "small-rna-plant-special.md")
    result.assert_file_exists(ref_path, "small-rna-plant-special.md exists")
    result.assert_file_not_empty(ref_path, "small-rna-plant-special.md is not empty")

    with open(ref_path, "r") as f:
        content = f.read()

    result.assert_true(len(content.split("\n")) >= 50, "Has at least 50 lines")

    # Plant-specific content
    plant_tokens = ["植物", "miRNA", "siRNA", "miRBase", "降解组", "拟南芥", "水稻"]
    for token in plant_tokens:
        result.assert_true(
            token in content,
            f"Contains plant-specific token: '{token}'"
        )

    # Key resources
    resources = ["psRNATarget", "miRBase"]
    for res in resources:
        result.assert_true(res in content, f"References resource: {res}")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
