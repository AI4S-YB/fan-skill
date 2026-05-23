#!/usr/bin/env python3
"""L2 Test: SKILL.md metadata validation."""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: SKILL.md Metadata")

    skill_path = os.path.join(SKILL_DIR, "SKILL.md")
    with open(skill_path, "r") as f:
        content = f.read()

    result.assert_true("skill_name: bio-plant-chipseq" in content, "Has skill_name")
    result.assert_true("version:" in content, "Has version")
    result.assert_true("category: plant-bioinformatics" in content, "Has category")
    result.assert_true("chip" in content.lower(), "References ChIP-seq")
    result.assert_true("## Overview" in content or "## When to Use" in content, "Has overview section")
    result.assert_true("## Decision Matrix" in content, "Has Decision Matrix section")
    result.assert_true("## Tool Catalog" in content, "Has Tool Catalog section")
    result.assert_true("## References" in content, "Has References section")
    result.assert_true("## Fallback" in content, "Has Fallback section")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
