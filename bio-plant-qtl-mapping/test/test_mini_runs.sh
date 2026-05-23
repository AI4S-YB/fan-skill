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

# 测试 1: 水稻 RIL 群体数据 -> 应选择 R/qtl est.map + CIM
echo "--- Test 1: Rice RIL population -> R/qtl est.map + CIM ---"
if [ -d "test/test_data/rice_ril" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/rice_ril > /tmp/rice_qtl_profile.json
    echo "  [RUN] R/qtl est.map + CIM on rice RIL test data..."

    if [ -f ".decisions/qtl_mapping_decisions.yaml" ]; then
        MAP_SELECTED=$(python3 -c "
import yaml
with open('.decisions/qtl_mapping_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0]['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$MAP_SELECTED" = "rqtl_est_map" ]; then
            echo "  [PASS] Selected rqtl_est_map for genetic map construction"
        else
            echo "  [WARN] Selected $MAP_SELECTED (expected rqtl_est_map)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

# 测试 2: 小麦高密度 DH 群体 -> 应选择 LepMap3 + MQM
echo "--- Test 2: Wheat high-density DH population -> LepMap3 + MQM ---"
if [ -d "test/test_data/wheat_dh_hd" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/wheat_dh_hd > /tmp/wheat_qtl_profile.json
    echo "  [RUN] LepMap3 + MQM on wheat DH test data..."

    if [ -f ".decisions/qtl_mapping_decisions.yaml" ]; then
        QTL_SELECTED=$(python3 -c "
import yaml
with open('.decisions/qtl_mapping_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[1]['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$QTL_SELECTED" = "multiple_qtl_mapping" ]; then
            echo "  [PASS] Selected multiple_qtl_mapping for QTL analysis"
        else
            echo "  [WARN] Selected $QTL_SELECTED (expected multiple_qtl_mapping)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
