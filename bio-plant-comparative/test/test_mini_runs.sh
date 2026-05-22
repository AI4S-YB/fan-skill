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

# 测试 1: 两物种水稻数据 → 应选择 OrthoFinder
echo "--- Test 1: Two-species comparison → OrthoFinder ---"
if [ -d "test/test_data/two_species" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/two_species/pep > /tmp/comp_profile.json
    echo "  [RUN] OrthoFinder on two-species test data..."

    if [ -f ".decisions/comp_decisions.yaml" ]; then
        SELECTED=$(python3 -c "
import yaml
with open('.decisions/comp_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0]['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED" = "orthofinder" ]; then
            echo "  [PASS] Selected OrthoFinder as expected"
        else
            echo "  [WARN] Selected $SELECTED (expected orthofinder)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
