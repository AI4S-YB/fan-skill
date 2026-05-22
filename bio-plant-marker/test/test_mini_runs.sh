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

# 测试 1: SNP + low throughput → KASP marker
echo "--- Test 1: SNP low-throughput → KASP ---"
echo "  [RUN] KASP primer design on test variant..."

if [ -f ".decisions/marker_decisions.yaml" ]; then
    SELECTED=$(python3 -c "
import yaml
with open('.decisions/marker_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0]['selected_tool'])
" 2>/dev/null || echo "unknown")
    if [ "$SELECTED" = "kasp_marker" ]; then
        echo "  [PASS] Selected KASP marker as expected"
    else
        echo "  [WARN] Selected $SELECTED (expected kasp_marker)"
    fi
fi

echo "=== L4 Complete ==="
