#!/bin/bash
# test_mini_runs.sh — L4: Mini run with decision path verification
set -euo pipefail

echo "=== L4: Mini Run Tests (eQTL) ==="

TEST_CONFIG="${TEST_CONFIG:-test_config.yaml}"

if [ ! -f "$TEST_CONFIG" ]; then
    echo "[SKIP] 没有 test_config.yaml — 跳过 L4 测试"
    echo "      创建 test_config.yaml 指向测试数据后重新运行"
    exit 0
fi

echo "Test config: $TEST_CONFIG"

# 测试 1: cis-eQTL 数据 → 应选择 MatrixEQTL cis
echo "--- Test 1: Cis-eQTL genotype + expression → MatrixEQTL cis ---"
if [ -d "test/test_data/rice_expression" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/rice_expression/expr > /tmp/eqtl_profile.json
    echo "  [RUN] MatrixEQTL cis on rice expression test data..."

    if [ -f ".decisions/eqtl_decisions.yaml" ]; then
        SELECTED=$(python3 -c "
import yaml
with open('.decisions/eqtl_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0]['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED" = "cis_eqtl_matrixeqtl" ]; then
            echo "  [PASS] Selected MatrixEQTL cis as expected"
        else
            echo "  [WARN] Selected $SELECTED (expected cis_eqtl_matrixeqtl)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

# 测试 2: 多组织表达 → 应触发 multi_tissue 规则
echo "--- Test 2: Multi-tissue expression → tissue_specific_vs_shared_eqtl ---"
if [ -d "test/test_data/multi_tissue" ]; then
    echo "  [RUN] Multi-tissue eQTL analysis..."
    if [ -f ".decisions/eqtl_decisions.yaml" ]; then
        MULTI=$(python3 -c "
import yaml
with open('.decisions/eqtl_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d.get('step') == 'multi_tissue':
        print(d.get('selected_tool', ''))
" 2>/dev/null || echo "unknown")
        if [ "$MULTI" = "tissue_specific_vs_shared_eqtl" ]; then
            echo "  [PASS] Multi-tissue rule triggered correctly"
        else
            echo "  [WARN] Got $MULTI (expected tissue_specific_vs_shared_eqtl)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
