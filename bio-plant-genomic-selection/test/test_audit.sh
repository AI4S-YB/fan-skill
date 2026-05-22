#!/bin/bash
set -euo pipefail
echo "=== L5: Decision Audit (GS) ==="
DECISIONS_DIR=".decisions"
if [ ! -d "$DECISIONS_DIR" ]; then
    echo "[SKIP] 无决策日志 — 先运行分析"
    exit 0
fi
REQUIRED_STEPS=("genotype_preprocessing" "model_selection" "cross_validation" "accuracy_metrics" "selection_strategy")
for step in "${REQUIRED_STEPS[@]}"; do
    if grep -q "$step" "$DECISIONS_DIR"/*.yaml 2>/dev/null; then
        echo "  [OK] Decision logged: $step"
    else
        echo "  [WARN] Missing decision: $step"
    fi
done
echo "=== L5 Complete ==="
