#!/bin/bash
# test_audit.sh — L5: Decision audit
set -euo pipefail

echo "=== L5: Decision Audit ==="

DECISIONS_DIR=".decisions"

if [ ! -d "$DECISIONS_DIR" ]; then
    echo "[SKIP] 无决策日志 — 先运行分析"
    exit 0
fi

REQUIRED_STEPS=("preprocessing" "clustering" "stage_specificity" "trend_analysis" "visualization")

for step in "${REQUIRED_STEPS[@]}"; do
    if grep -q "$step" "$DECISIONS_DIR"/*.yaml 2>/dev/null; then
        echo "  [OK] Decision logged: $step"
    else
        echo "  [WARN] Missing decision: $step"
    fi
done

# 检查 rule 模式决策有 rule_id
python3 << 'PYEOF'
import yaml, os, sys

decisions_dir = ".decisions"
issues = 0

for f in os.listdir(decisions_dir):
    if not f.endswith('.yaml'):
        continue
    with open(os.path.join(decisions_dir, f)) as fh:
        data = yaml.safe_load(fh)
    for d in data.get('decisions', []):
        if d['mode'] == 'rule' and 'rule_id' not in d:
            print(f"  [ISSUE] Rule decision without rule_id: {d['step']}")
            issues += 1
        if d['mode'] in ('expert', 'hybrid_fallback') and 'reasoning' not in d:
            print(f"  [ISSUE] Expert decision without reasoning: {d['step']}")
            issues += 1

if issues == 0:
    print("  [OK] All decisions properly documented")
else:
    print(f"  [WARN] {issues} decisions with incomplete documentation")

sys.exit(0)
PYEOF

echo "=== L5 Complete ==="
