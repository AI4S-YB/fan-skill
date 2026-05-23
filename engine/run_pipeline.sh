#!/bin/bash
# run_pipeline.sh — marker-based checkpoint pipeline for long-running analyses
set -euo pipefail

CHECKPOINT_DIR=".checkpoints"
PROGRESS_FILE="PROGRESS.md"
PID_FILE=".pipeline.pid"

marker_file() { echo "$CHECKPOINT_DIR/${1}.done"; }

run_step() {
    local step=$1; shift
    local marker
    marker=$(marker_file "$step")
    mkdir -p "$CHECKPOINT_DIR"
    if [ -f "$marker" ]; then
        echo "[SKIP] $step — completed at $(cat "$marker")"
        return 0
    fi
    echo "[RUN]  $step — $(date '+%H:%M:%S')"
    if "$@"; then
        date '+%Y-%m-%d %H:%M:%S' > "$marker"
        echo "[DONE] $step"
    else
        echo "[FAIL] $step — delete $marker to retry"
        return 1
    fi
}

write_progress() {
    local msg="$1"
    echo "$(date '+%H:%M:%S') — $msg" >> "$PROGRESS_FILE"
}

# PID lock
if [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
    echo "Pipeline already running (PID=$(cat $PID_FILE))"
    exit 1
fi
echo $$ > "$PID_FILE"
trap "rm -f $PID_FILE" EXIT

echo "# Pipeline Progress — $(date)" > "$PROGRESS_FILE"
echo "Usage: delete .checkpoints/<step>.done to re-run that step" >> "$PROGRESS_FILE"
