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

    ref_path = os.path.join(SKILL_DIR, "references", "chipseq-plant-special.md")
    result.assert_file_exists(ref_path, "chipseq-plant-special.md exists")
    result.assert_file_not_empty(ref_path, "chipseq-plant-special.md is not empty")

    with open(ref_path, "r") as f:
        content = f.read()

    result.assert_true(len(content.split("\n")) >= 50, "Has at least 50 lines")

    plant_tokens = ["植物", "组蛋白", "基因组", "拟南芥", "抗体", "转录因子", "重复序列"]
    for token in plant_tokens:
        result.assert_true(token in content, f"Contains plant-specific token: '{token}'")

    resources = ["MACS2", "DiffBind", "ChIPseeker"]
    for res in resources:
        result.assert_true(res in content, f"References resource: {res}")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
