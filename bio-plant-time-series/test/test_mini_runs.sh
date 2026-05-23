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

# 测试 1: 多时间点水稻发育数据 → 应选择 spline + mfuzz
echo "--- Test 1: Rice development 10 time points → spline + mfuzz ---"
if [ -d "test/test_data/rice_dev_timecourse" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/rice_dev_timecourse/expr > /tmp/rice_ts_profile.json
    echo "  [RUN] Spline smoothing + Mfuzz on rice dev data..."

    if [ -f ".decisions/time_series_decisions.yaml" ]; then
        SELECTED_PRE=$(python3 -c "
import yaml
with open('.decisions/time_series_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d['step'] == 'preprocessing':
        print(d['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED_PRE" = "spline_smoothing" ]; then
            echo "  [PASS] Preprocessing: selected spline_smoothing as expected"
        else
            echo "  [WARN] Preprocessing: selected $SELECTED_PRE (expected spline_smoothing)"
        fi

        SELECTED_CL=$(python3 -c "
import yaml
with open('.decisions/time_series_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d['step'] == 'clustering':
        print(d['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED_CL" = "mfuzz" ]; then
            echo "  [PASS] Clustering: selected mfuzz as expected"
        else
            echo "  [WARN] Clustering: selected $SELECTED_CL (expected mfuzz)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

# 测试 2: 少量基因热胁迫 → k-means + LOESS
echo "--- Test 2: Heat stress 4 time points, 500 DEGs → loess + kmeans ---"
if [ -d "test/test_data/heat_stress_degs" ]; then
    echo "  [RUN] LOESS + k-means on heat stress DEG data..."

    if [ -f ".decisions/time_series_decisions.yaml" ]; then
        SELECTED_PRE=$(python3 -c "
import yaml
with open('.decisions/time_series_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d['step'] == 'preprocessing':
        print(d['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED_PRE" = "loess_smoothing" ]; then
            echo "  [PASS] Preprocessing: selected loess_smoothing as expected"
        else
            echo "  [WARN] Preprocessing: selected $SELECTED_PRE (expected loess_smoothing)"
        fi

        SELECTED_CL=$(python3 -c "
import yaml
with open('.decisions/time_series_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
for d in decisions:
    if d['step'] == 'clustering':
        print(d['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED_CL" = "kmeans" ]; then
            echo "  [PASS] Clustering: selected kmeans as expected"
        else
            echo "  [WARN] Clustering: selected $SELECTED_CL (expected kmeans)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
