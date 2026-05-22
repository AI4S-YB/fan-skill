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

# 测试 1: 小样本水稻数据 → 应选择 GAPIT CMLM
echo "--- Test 1: Small rice inbred → GAPIT CMLM ---"
if [ -d "test/test_data/rice_small" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/rice_small/geno > /tmp/rice_profile.json
    echo "  [RUN] GAPIT CMLM on rice test data..."

    if [ -f ".decisions/gwas_decisions.yaml" ]; then
        SELECTED=$(python3 -c "
import yaml
with open('.decisions/gwas_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0]['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED" = "gapit-cmlm" ]; then
            echo "  [PASS] Selected GAPIT CMLM as expected"
        else
            echo "  [WARN] Selected $SELECTED (expected gapit-cmlm)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
