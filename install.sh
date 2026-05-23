#!/bin/bash
# install.sh — install fan-skill v2.0 to Claude Code
set -euo pipefail

SKILLS_DIR="${CLAUDE_SKILLS_DIR:-$HOME/.claude/skills}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INSTALL_DIR="$SKILLS_DIR/fan-skill"

echo "=== Fan-Skill v2.0 Install ==="
echo "Install target: $INSTALL_DIR"

mkdir -p "$INSTALL_DIR"

# Copy v2.0 structure
for component in SKILL.md engine knowledge-base tool-registry references templates theme; do
    if [ -e "$SCRIPT_DIR/$component" ]; then
        echo "  Installing: $component"
        cp -rp "$SCRIPT_DIR/$component" "$INSTALL_DIR/$component"
    fi
done

echo "=== Installation complete ==="
echo "fan-skill installed to: $INSTALL_DIR"
echo ""
echo "To use: just describe your biological question in Claude Code."
echo "fan-skill will automatically match your intent to the knowledge base."
