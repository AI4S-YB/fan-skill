#!/bin/bash
# install-opencode.sh — install fan-skill for OpenCode (project-level only)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INSTALL_DIR="$SCRIPT_DIR/.opencode/skill/fan-skill"

echo "=== Fan-Skill Install (OpenCode) ==="
echo "Install target: $INSTALL_DIR"

rm -rf "$INSTALL_DIR"
mkdir -p "$(dirname "$INSTALL_DIR")"

for component in SKILL.md engine knowledge-base tool-registry references templates theme; do
    if [ -e "$SCRIPT_DIR/$component" ]; then
        echo "  Installing: $component"
        cp -rp "$SCRIPT_DIR/$component" "$INSTALL_DIR/$component"
    fi
done

echo "=== Installation complete ==="
echo "fan-skill installed to: $INSTALL_DIR"
echo "To use: just describe your biological question in OpenCode."
