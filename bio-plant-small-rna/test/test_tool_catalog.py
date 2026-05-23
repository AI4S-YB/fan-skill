#!/usr/bin/env python3
"""L2 Test: Tool catalog completeness and content."""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CATALOG_DIR = os.path.join(SKILL_DIR, "tool-catalog")

def main():
    result = TestResult()
    print("Test: Tool Catalog Completeness")

    expected_files = [
        "mirna-prediction.md",
        "target-prediction.md",
        "degradome.md",
        "differential-mirna.md",
        "preprocessing.md",
        "visualization.md",
    ]

    content_checks = {
        "mirna-prediction.md": ["miRDeep2", "miRNA", "miRBase"],
        "target-prediction.md": ["psRNATarget", "靶基因", "miRNA"],
        "degradome.md": ["CleaveLand", "降解组", "PARE"],
        "differential-mirna.md": ["DESeq2", "差异", "edgeR"],
        "preprocessing.md": ["FastQC", "Cutadapt", "18", "30"],
        "visualization.md": ["Volcano", "Heatmap", "ggplot2"],
    }

    for f in expected_files:
        path = os.path.join(CATALOG_DIR, f)
        result.assert_file_exists(path, f"File '{f}' exists")
        result.assert_file_not_empty(path, f"File '{f}' is not empty")

        with open(path, "r") as fh:
            content = fh.read()

        result.assert_true(len(content.split("\n")) >= 10, f"'{f}' has >= 10 lines")
        result.assert_true(content.startswith("#"), f"'{f}' starts with heading")

        if f in content_checks:
            for keyword in content_checks[f]:
                result.assert_true(
                    keyword.lower() in content.lower(),
                    f"'{f}' contains keyword '{keyword}'"
                )

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
