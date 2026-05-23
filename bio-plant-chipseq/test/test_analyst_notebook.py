#!/usr/bin/env python3
"""L2 Test: Analyst notebook content validation."""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Analyst Notebook")

    nb_path = os.path.join(SKILL_DIR, "analyst-notebook.md")
    result.assert_file_exists(nb_path, "analyst-notebook.md exists")
    result.assert_file_not_empty(nb_path, "analyst-notebook.md is not empty")

    with open(nb_path, "r") as f:
        content = f.read()

    result.assert_true(len(content.split("\n")) >= 80, "Has at least 80 lines")
    result.assert_true("植物" in content, "Contains Chinese content")

    pipeline_stages = [
        ("FastQC", "Quality control"),
        ("Bowtie2", "Read alignment"),
        ("MACS2", "Peak calling"),
        ("DiffBind", "Differential binding"),
        ("ChIPseeker", "Peak annotation"),
        ("MEME", "Motif analysis"),
    ]
    for token, label in pipeline_stages:
        result.assert_true(token in content, f"Covers {label} ({token})")

    result.assert_true("```" in content, "Contains code blocks")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
