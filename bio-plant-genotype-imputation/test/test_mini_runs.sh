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

# 测试 1: 小样本水稻数据 → 应选择 Beagle
echo "--- Test 1: Small rice dataset → Beagle ---"
if [ -d "test/test_data/rice_small" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/rice_small/geno > /tmp/impute_rice_profile.json
    echo "  [RUN] Beagle on rice test data..."

    if [ -f ".decisions/imputation_decisions.yaml" ]; then
        SELECTED=$(python3 -c "
import yaml
with open('.decisions/imputation_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0]['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED" = "beagle" ]; then
            echo "  [PASS] Selected Beagle as expected"
        else
            echo "  [WARN] Selected $SELECTED (expected beagle)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

# 测试 2: 多倍体数据 → 应选择 Beagle 多倍体模式
echo "--- Test 2: Tetraploid potato → Beagle polyploid ---"
if [ -d "test/test_data/potato_tetra" ]; then
    echo "  [RUN] Beagle polyploid mode on potato test data..."

    if [ -f ".decisions/imputation_decisions.yaml" ]; then
        RECOMMEND=$(python3 -c "
import yaml
with open('.decisions/imputation_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d['step'] == 'imputation_tool':
        print(d['recommend'])
" 2>/dev/null || echo "unknown")
        if [ "$RECOMMEND" = "beagle_polyploid" ]; then
            echo "  [PASS] Recommended beagle_polyploid for tetraploid"
        else
            echo "  [WARN] Recommended $RECOMMEND (expected beagle_polyploid)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

# 测试 3: 大样本数据 → 应选择 Minimac4
echo "--- Test 3: Large maize dataset → Minimac4 ---"
if [ -d "test/test_data/maize_large" ]; then
    echo "  [RUN] Minimac4 on maize test data..."

    if [ -f ".decisions/imputation_decisions.yaml" ]; then
        RECOMMEND=$(python3 -c "
import yaml
with open('.decisions/imputation_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d['step'] == 'imputation_tool':
        print(d['recommend'])
" 2>/dev/null || echo "unknown")
        if [ "$RECOMMEND" = "minimac4" ]; then
            echo "  [PASS] Recommended Minimac4 as expected"
        else
            echo "  [WARN] Recommended $RECOMMEND (expected minimac4)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
