#!/bin/bash
# install-gemini.sh — install fan-skill for Gemini CLI
set -euo pipefail

SKILLS_DIR="${GEMINI_SKILLS_DIR:-$HOME/.gemini/skills}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INSTALL_DIR="$SKILLS_DIR/fan-skill"

echo "=== Fan-Skill Install (Gemini CLI) ==="
echo "Install target: $INSTALL_DIR"

rm -rf "$INSTALL_DIR"
mkdir -p "$INSTALL_DIR"

for component in SKILL.md engine knowledge-base tool-registry references templates theme; do
    if [ -e "$SCRIPT_DIR/$component" ]; then
        echo "  Installing: $component"
        cp -rp "$SCRIPT_DIR/$component" "$INSTALL_DIR/$component"
    fi
done

echo "=== Installation complete ==="
echo "fan-skill installed to: $INSTALL_DIR"
echo "To use: just describe your biological question in Gemini CLI."
