#!/bin/bash
# log_decision.sh — record analysis decisions to log
# Usage: bash log_decision.sh --step <name> --mode <rule|expert|hybrid_fallback> \
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

# Infer Skill name from output_dir
SKILL_NAME=$(basename "$(dirname "$OUTPUT_DIR")")
LOG_FILE="${OUTPUT_DIR}/${SKILL_NAME}_decisions.yaml"

# Initialize log file
if [ ! -f "$LOG_FILE" ]; then
    echo "decisions:" > "$LOG_FILE"
fi

TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

# Append decision record
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
    # Python yaml not available — fallback
    echo "  - step: $STEP" >> "$LOG_FILE"
    echo "    mode: $MODE" >> "$LOG_FILE"
    echo "    timestamp: $TIMESTAMP" >> "$LOG_FILE"
    echo "    selected_tool: $SELECTED" >> "$LOG_FILE"
    echo "    reasoning: \"$REASON\"" >> "$LOG_FILE"
}

echo "[LOG] $STEP → $SELECTED ($MODE)"
