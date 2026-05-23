#!/usr/bin/env python3
"""L2 Test: Structure validation for bio-plant-small-rna."""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Structure Validation")

    required_dirs = ["tool-catalog", "references", "test"]
    for d in required_dirs:
        result.assert_file_exists(os.path.join(SKILL_DIR, d), f"Directory '{d}' exists")

    required_files = ["SKILL.md", "decision-matrix.yaml", "analyst-notebook.md"]
    for f in required_files:
        path = os.path.join(SKILL_DIR, f)
        result.assert_file_exists(path, f"File '{f}' exists")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
