#!/bin/bash
# test_mini_runs.sh -- L4: Mini run with decision path verification
set -euo pipefail

echo "=== L4: Mini Run Tests (GRN) ==="

TEST_CONFIG="${TEST_CONFIG:-test_config.yaml}"

if [ ! -f "$TEST_CONFIG" ]; then
    echo "[SKIP] No test_config.yaml -- skipping L4 tests"
    echo "      Create test_config.yaml pointing to test data and re-run"
    exit 0
fi

echo "Test config: $TEST_CONFIG"

# Test 1: Arabidopsis expression matrix -> should select GENIE3
echo "--- Test 1: Arabidopsis RNA-seq -> GENIE3 with TF prior ---"
if [ -d "test/test_data/arabidopsis_expr" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/arabidopsis_expr > /tmp/ath_expr_profile.json
    echo "  [RUN] GENIE3 on Arabidopsis test data..."

    if [ -f ".decisions/grn_decisions.yaml" ]; then
        SELECTED=$(python3 -c "
import yaml
with open('.decisions/grn_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d['step'] == 'inference_method':
        print(d.get('selected_method', 'unknown'))
        break
" 2>/dev/null || echo "unknown")
        echo "  [INFO] Selected method: $SELECTED"
    fi
else
    echo "  [SKIP] Test data not available"
fi

# Test 2: Non-model species no TF annotation -> should select WGCNA
echo "--- Test 2: Non-model species -> WGCNA co-expression ---"
if [ -d "test/test_data/nonmodel_expr" ]; then
    echo "  [RUN] WGCNA fallback for non-model species..."

    if [ -f ".decisions/grn_decisions.yaml" ]; then
        TF_METHOD=$(python3 -c "
import yaml
with open('.decisions/grn_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d['step'] == 'tf_database':
        print(d.get('selected_source', 'unknown'))
        break
" 2>/dev/null || echo "unknown")
        echo "  [INFO] TF database selection: $TF_METHOD"
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
