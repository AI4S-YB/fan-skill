#!/bin/bash
# qc_check.sh — 根据 qc-thresholds.yaml 检查 QC 指标
# 用法: bash qc_check.sh --analysis <type> --step <name> --metric <name> --value <number>

ANALYSIS=""
STEP=""
METRIC=""
VALUE=""
THRESHOLDS_FILE="$(dirname "$0")/../references/qc-thresholds.yaml"
QC_REPORT=".qc_report.txt"

while [[ $# -gt 0 ]]; do
    case $1 in
        --analysis) ANALYSIS="$2"; shift 2 ;;
        --step) STEP="$2"; shift 2 ;;
        --metric) METRIC="$2"; shift 2 ;;
        --value) VALUE="$2"; shift 2 ;;
        *) shift ;;
    esac
done

if [ -z "$ANALYSIS" ] || [ -z "$METRIC" ] || [ -z "$VALUE" ]; then
    echo "Usage: qc_check.sh --analysis <type> --step <name> --metric <name> --value <number>"
    exit 1
fi

# 从 YAML 中提取阈值
THRESHOLD=$(grep -A3 "$METRIC:" "$THRESHOLDS_FILE" 2>/dev/null | grep "threshold:" | grep -oP '[0-9.]+' | head -1)

if [ -z "$THRESHOLD" ]; then
    echo "[INFO] $STEP/$METRIC=$VALUE — 无预设阈值，跳过自动判定"
    echo "$(date): $ANALYSIS/$STEP/$METRIC=$VALUE INFO (no threshold)" >> "$QC_REPORT"
    exit 0
fi

# 浮点比较
RESULT=$(python3 -c "
value = float('$VALUE')
threshold = float('$THRESHOLD')
if value >= threshold:
    print('PASS')
elif value >= threshold * 0.8:
    print('WARN')
else:
    print('FAIL')
" 2>/dev/null || echo "UNKNOWN")

echo "[$RESULT] $STEP/$METRIC=$VALUE (threshold=$THRESHOLD)"
echo "$(date): $ANALYSIS/$STEP/$METRIC=$VALUE $RESULT (threshold=$THRESHOLD)" >> "$QC_REPORT"
