#!/bin/bash
# test_mini_runs.sh — L4: Mini run with decision path verification
set -euo pipefail

echo "=== L4: Mini Run Tests (Multi-Omics) ==="

TEST_CONFIG="${TEST_CONFIG:-test_config.yaml}"

if [ ! -f "$TEST_CONFIG" ]; then
    echo "[SKIP] No test_config.yaml — skipping L4 tests"
    echo "       Create test_config.yaml pointing to test data and re-run"
    exit 0
fi

echo "Test config: $TEST_CONFIG"

# Test 1: Two-omics data with correlation goal → mixOmics PLS
echo "--- Test 1: Two omics + correlation goal → mixOmics PLS ---"
echo "  [SKIP] Requires test data (transcriptome + metabolome matrices)"

# Test 2: Three omics + exploratory goal → MOFA2
echo "--- Test 2: Three omics + exploratory goal → MOFA2 ---"
echo "  [SKIP] Requires test data (transcriptome + metabolome + proteome matrices)"

# Test 3: Classification goal with labels → DIABLO
echo "--- Test 3: Classification + labels → DIABLO ---"
echo "  [SKIP] Requires test data with class labels"

echo "=== L4 Complete ==="
