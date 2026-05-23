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
        "peak-calling.md",
        "diff-binding.md",
        "motif-enrichment.md",
        "annotation.md",
        "preprocessing.md",
        "visualization.md",
    ]

    content_checks = {
        "peak-calling.md": ["MACS2", "narrow", "broad", "转录因子", "组蛋白"],
        "diff-binding.md": ["DiffBind", "DESeq2", "差异", "edgeR"],
        "motif-enrichment.md": ["MEME", "motif", "JASPAR", "HOMER"],
        "annotation.md": ["ChIPseeker", "peak", "基因", "promoter"],
        "preprocessing.md": ["FastQC", "Bowtie2", "SAMtools", "植物"],
        "visualization.md": ["deepTools", "IGV", "bigWig", "Heatmap"],
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
            found_count = 0
            for keyword in content_checks[f]:
                if keyword.lower() in content.lower():
                    found_count += 1
            result.assert_true(found_count >= 3, f"'{f}' contains at least 3 key terms (found {found_count})")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
