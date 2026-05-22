#!/bin/bash
# test_audit.sh — L5: Decision audit
set -euo pipefail

echo "=== L5: Decision Audit ==="

DECISIONS_DIR=".decisions"

if [ ! -d "$DECISIONS_DIR" ]; then
    echo "[SKIP] 无决策日志 — 先运行分析"
    exit 0
fi

REQUIRED_STEPS=("ld_pruning" "pca" "admixture" "fst" "phylogeny")

for step in "${REQUIRED_STEPS[@]}"; do
    if grep -q "$step" "$DECISIONS_DIR"/*.yaml 2>/dev/null; then
        echo "  [OK] Decision logged: $step"
    else
        echo "  [WARN] Missing decision: $step"
    fi
done

echo "=== L5 Complete ==="
