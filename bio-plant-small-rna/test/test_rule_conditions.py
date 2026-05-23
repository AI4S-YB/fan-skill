#!/usr/bin/env python3
"""L2 Test: Decision matrix - species detection conditions."""
import os
import sys
import yaml
sys.path.insert(0, os.path.dirname(__file__))
from test_engine import TestResult

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    result = TestResult()
    print("Test: Species Detection Conditions")

    matrix_path = os.path.join(SKILL_DIR, "decision-matrix.yaml")
    with open(matrix_path, "r") as f:
        dm = yaml.safe_load(f)

    mirna_rules = dm["mirna_prediction"]

    # Rule mirna-plant-species-001: model species
    r001 = [r for r in mirna_rules if r["rule_id"] == "mirna-plant-species-001"][0]
    species = r001["condition"]["species_detected"]

    model_species = ["oryza_sativa", "zea_mays", "arabidopsis_thaliana", "glycine_max"]
    result.assert_equal(len(species), 4, "Exactly 4 model species")
    for sp in model_species:
        result.assert_true(sp in species, f"Model species list includes {sp}")

    # Rule mirna-novel-002: unknown species
    r002 = [r for r in mirna_rules if r["rule_id"] == "mirna-novel-002"][0]
    result.assert_equal(r002["condition"]["species_detected"], "unknown", "Novel rule targets unknown species")

    # Priority ordering: model species > novel
    result.assert_true(r001["priority"] > r002["priority"],
                       f"Model species priority ({r001['priority']}) > novel priority ({r002['priority']})")

    return result.summary()

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
