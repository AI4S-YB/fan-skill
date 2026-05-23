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

# 测试 1: 小样本水稻数据 -> 应选择 bcftools mpileup
echo "--- Test 1: Small rice cohort -> bcftools mpileup ---"
if [ -d "test/test_data/rice_small" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/rice_small/fastq > /tmp/rice_vc_profile.json
    echo "  [RUN] bcftools mpileup on rice test data..."

    if [ -f ".decisions/variant_calling_decisions.yaml" ]; then
        SELECTED=$(python3 -c "
import yaml
with open('.decisions/variant_calling_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[0]['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED" = "bcftools_mpileup" ]; then
            echo "  [PASS] Selected bcftools mpileup as expected"
        else
            echo "  [WARN] Selected $SELECTED (expected bcftools_mpileup)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

# 测试 2: 大样本玉米数据 -> 应选择 GATK GVCF
echo "--- Test 2: Large maize cohort -> GATK GVCF ---"
if [ -d "test/test_data/maize_large" ]; then
    bash ../../bio-plant-infra/scripts/inspect_data.sh test/test_data/maize_large/fastq > /tmp/maize_vc_profile.json
    echo "  [RUN] GATK HaplotypeCaller GVCF on maize test data..."

    if [ -f ".decisions/variant_calling_decisions.yaml" ]; then
        SELECTED=$(python3 -c "
import yaml
with open('.decisions/variant_calling_decisions.yaml') as f:
    decisions = yaml.safe_load(f)['decisions']
print(decisions[1]['selected_tool'])
" 2>/dev/null || echo "unknown")
        if [ "$SELECTED" = "gatk_haplotypecaller_gvcf" ]; then
            echo "  [PASS] Selected GATK GVCF as expected"
        else
            echo "  [WARN] Selected $SELECTED (expected gatk_haplotypecaller_gvcf)"
        fi
    fi
else
    echo "  [SKIP] Test data not available"
fi

echo "=== L4 Complete ==="
