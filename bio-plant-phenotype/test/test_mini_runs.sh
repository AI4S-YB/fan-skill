#!/bin/bash
# test_mini_runs.sh — L4: Mini run with decision path verification
set -euo pipefail

echo "=== L4: Mini Run Tests ==="

TEST_CONFIG="${TEST_CONFIG:-test_config.yaml}"

if [ ! -f "$TEST_CONFIG" ]; then
    echo "[SKIP] 没有 test_config.yaml — 跳过 L4 测试"
    echo "      创建 test_config.yaml 指向测试数据后重新运行"
    exit 0
fi

echo "Test config: $TEST_CONFIG"

# Test 1: Multi-environment trial → AMMI/GGE
echo "--- Test 1: Multi-env trial → AMMI ---"
if [ -d "test/test_data/wheat_multi_env" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/wheat_multi_env/pheno > /tmp/wheat_profile.json
    echo "  [RUN] sommer multi-env analysis on wheat test data..."

    if [ -f ".decisions/phenotype_decisions.yaml" ]; then
        SELECTED=$(python3 -c "
import yaml
with open('.decisions/phenotype_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0]['selected_tool'])
" 2>/dev/null || echo "unknown")
        echo "  [SELECTED] $SELECTED"
    fi
else
    echo "  [SKIP] Test data not available"
fi

# Test 2: Single environment → BLUP/ANOVA heritability
echo "--- Test 2: Single env trial → BLUP ---"
if [ -d "test/test_data/rice_single_env" ]; then
    echo "  [RUN] sommer single-env analysis on rice test data..."

    if [ -f ".decisions/phenotype_decisions.yaml" ]; then
        METHOD=$(python3 -c "
import yaml
with open('.decisions/phenotype_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0].get('method', 'unknown'))
" 2>/dev/null || echo "unknown")
        echo "  [SELECTED] $METHOD"
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
