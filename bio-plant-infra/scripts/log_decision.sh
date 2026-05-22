#!/bin/bash
# log_decision.sh — 记录分析决策到日志
# 用法: bash log_decision.sh --step <name> --mode <rule|expert|hybrid_fallback> \
#           --selected <tool> --reason "<text>" [--overridden_from <rule_id>] \
#           [--output_dir <path>]

OUTPUT_DIR=".decisions"
STEP=""
MODE=""
SELECTED=""
REASON=""
OVERRIDDEN=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --step) STEP="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        --selected) SELECTED="$2"; shift 2 ;;
        --reason) REASON="$2"; shift 2 ;;
        --overridden_from) OVERRIDDEN="$2"; shift 2 ;;
        --output_dir) OUTPUT_DIR="$2/.decisions"; shift 2 ;;
        *) shift ;;
    esac
done

if [ -z "$STEP" ] || [ -z "$MODE" ]; then
    echo "Usage: log_decision.sh --step <name> --mode <rule|expert|hybrid_fallback> --selected <tool> --reason <text>"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# 从 output_dir 推断 Skill 名
SKILL_NAME=$(basename "$(dirname "$OUTPUT_DIR")")
LOG_FILE="${OUTPUT_DIR}/${SKILL_NAME}_decisions.yaml"

# 初始化日志文件
if [ ! -f "$LOG_FILE" ]; then
    echo "decisions:" > "$LOG_FILE"
fi

TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

# 追加决策记录
python3 -c "
import yaml, sys

entry = {
    'step': '$STEP',
    'mode': '$MODE',
    'timestamp': '$TIMESTAMP',
    'selected_tool': '$SELECTED',
    'reasoning': '$REASON',
}

if '$OVERRIDDEN':
    entry['rule_id_overridden'] = '$OVERRIDDEN'

with open('$LOG_FILE', 'a') as f:
    yaml.dump([entry], f, default_flow_style=False, allow_unicode=True)
" 2>/dev/null || {
    # Python yaml 不可用时的 fallback
    echo "  - step: $STEP" >> "$LOG_FILE"
    echo "    mode: $MODE" >> "$LOG_FILE"
    echo "    timestamp: $TIMESTAMP" >> "$LOG_FILE"
    echo "    selected_tool: $SELECTED" >> "$LOG_FILE"
    echo "    reasoning: \"$REASON\"" >> "$LOG_FILE"
}

echo "[LOG] $STEP → $SELECTED ($MODE)"
